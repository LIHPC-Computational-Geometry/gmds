#include <gmds/utils/OrientedGraph.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    GraphNode::GraphNode(const TCellID ANb)
    :m_id(ANb)
    {}
    /*----------------------------------------------------------------------------*/
    bool GraphNode::addOutEdge(GraphEdge* AEdge)
    {
        if(AEdge->tail()!=this)
            return false;
        
        m_out_edges.push_back(AEdge);
        return true;
    }
    /*----------------------------------------------------------------------------*/
    bool GraphNode::addInEdge(GraphEdge* AEdge)
    {
        if(AEdge->head()!=this)
            return false;
        
        m_in_edges.push_back(AEdge);
        return  true;
    }
    /*----------------------------------------------------------------------------*/
    TCellID GraphNode::id() const
    {
        return m_id;
    }
    
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*>& GraphNode::outEdges()
    {
        return m_out_edges;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*>& GraphNode::inEdges()
    {
        return m_in_edges;
    }
    /*----------------------------------------------------------------------------*/
    GraphEdge::GraphEdge(const TCellID AID,  GraphNode* ATail,  GraphNode* AHead)
    :m_id(AID), m_tail(ATail),m_head(AHead){
    }
    /*----------------------------------------------------------------------------*/
    TCellID GraphEdge::id()
    {
        return m_id;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*> GraphEdge::getEdgesStartingFrom(GraphNode* ANode){
        std::vector<GraphEdge*> edges;
        
        if(ANode!=m_head && ANode!=m_tail){
            return edges;
        }
        
        std::vector<GraphEdge*> all_edges = ANode->outEdges();
        for(auto e:all_edges){
            if(e->id()!=this->id()){
//                std::cout<<"Sharing edge: "<<e->tail()->id()<<", "
//                <<e->head()->id()<<std::endl;
                edges.push_back(e);
            }
        }
        
        return edges;
    }
    
    /*----------------------------------------------------------------------------*/
    GraphNode* GraphEdge::tail()
    {
        return m_tail;
    }
    /*----------------------------------------------------------------------------*/
    GraphNode* GraphEdge::head()
    {
        return m_head;
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::OrientedGraph(const int ANbNodes)
    {
        m_nodes.resize(ANbNodes);
        for(auto i=0; i<ANbNodes;i++){
            m_nodes[i].m_id=i;
            m_map_nodes[i]=&(m_nodes[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::OrientedGraph(const std::set<TCellID>& AIDS)
    {
        m_nodes.resize(AIDS.size());
        auto i=0;
        for(auto id:AIDS){
            m_nodes[i].m_id=id;
            m_map_nodes[id]=&(m_nodes[i]);
            i++;
        }
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::~OrientedGraph()
    {}
    /*----------------------------------------------------------------------------*/
    bool OrientedGraph::addEdge(const TCellID AID,
                                const TCellID ATail,
                                const TCellID AHead)
    {
        m_edges.push_back(GraphEdge(AID,
                                    m_map_nodes[ATail],
                                    m_map_nodes[AHead]));
        return true;
    }
    
    /*----------------------------------------------------------------------------*/
    void OrientedGraph::updateNodes()
    {
        
        for(int i=0; i<m_edges.size(); i++){
            GraphEdge* ei = &m_edges[i];
            ei->tail()->addOutEdge(ei);
            ei->head()->addInEdge(ei);
            
        }
    }
    /*----------------------------------------------------------------------------*/
    GraphNode* OrientedGraph::node(const TCellID AIndex)
    {
        if(m_map_nodes.find(AIndex)==m_map_nodes.end())
            throw GMDSException("Not found");
        
        return m_map_nodes[AIndex];
        
    }
    /*----------------------------------------------------------------------------*/
    GraphEdge* OrientedGraph::edge(const TCellID AIndex)
    {
        return &(m_edges[AIndex]);
    }
    /*----------------------------------------------------------------------------*/
    int OrientedGraph::nbNodes() const
    {
        return m_nodes.size();
    }
    /*----------------------------------------------------------------------------*/
    int OrientedGraph::nbEdges() const
    {
        return m_edges.size();
    }
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr,  gmds::OrientedGraph& AG)
{
    AStr<<"Oriented Graph ("<<AG.nbNodes()<<", "<<AG.nbEdges()<<")"<<std::endl;
    std::cout<<" Nodes: ";
    for(auto i:AG.m_nodes)
        std::cout<<i.id()<<" ";
    std::cout<<"\n Edges: ";
    for(auto i:AG.m_edges)
        std::cout<<"["<<i.head()->id()<<", "<<i.tail()->id()<<"] ";
    std::cout<<std::endl;
    return AStr;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
