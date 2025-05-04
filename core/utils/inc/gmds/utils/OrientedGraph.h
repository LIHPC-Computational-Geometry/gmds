/*----------------------------------------------------------------------------*/
#ifndef GMDS_ORIENTED_GRAPH_H_
#define GMDS_ORIENTED_GRAPH_H_
/*----------------------------------------------------------------------------*/
// STL Header files
#include <vector>
#include <set>
#include <map>
#include <iostream>
/*----------------------------------------------------------------------------*/
// GMDS Header files
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{

    //forward declaration
    class GraphEdge;
    class OrientedGraph;
    /*------------------------------------------------------------------------*/
    /** \class GraphNode
     *
     *  \brief Class defining a graph node. A graph node store a reference on 
     *         all the edges starting from it
     */
    /*------------------------------------------------------------------------*/
    class GraphNode{
        friend class OrientedGraph;
    public:
        /*--------------------------------------------------------------------*/
        /** \brief Default constructor
         */
        /*--------------------------------------------------------------------*/
        GraphNode(const TCellID ANb=NullID);

        /*--------------------------------------------------------------------*/
        /** \brief Add an edge starting from *this. Do nothing if \p AEdge does
         *         start from *this, i.e. AEdge.tail != *this
         *
         * \param[in] AEdge the edge to add.
         *
         * \return true if the edge was added, false otherwise
         */
        /*--------------------------------------------------------------------*/
        bool addOutEdge( GraphEdge* AEdge);
        
        /*--------------------------------------------------------------------*/
        /** \brief Add an edge arriving to *this. Do nothing if \p AEdge does
         *         end on *this, i.e. AEdge.head != *this
         *
         * \param[in] AEdge the edge to add.
         *
         * \return true if the edge was added, false otherwise
         */
        /*--------------------------------------------------------------------*/
        bool addInEdge( GraphEdge* AEdge);
        
        /*--------------------------------------------------------------------*/
        /** \brief Returns a reference onto the edges starting from *this
         */
        /*--------------------------------------------------------------------*/
        std::vector<GraphEdge*>& outEdges();
        /*--------------------------------------------------------------------*/
        /** \brief Returns a reference onto the edges coming to *this
         */
        /*--------------------------------------------------------------------*/
        std::vector<GraphEdge*>& inEdges();
        /*--------------------------------------------------------------------*/
        /** \brief Returns the node id
         */
        /*--------------------------------------------------------------------*/
        TCellID id() const;
    private:
        
        TCellID m_id;
        std::vector<GraphEdge*> m_out_edges;
        std::vector<GraphEdge*> m_in_edges;
    };
    
    /*------------------------------------------------------------------------*/
    /** \class GraphEdge
     *
     *  \brief Class defining a graph edge oriented from tail to head.
     */
    /*------------------------------------------------------------------------*/
    class GraphEdge{
    public:
        /*--------------------------------------------------------------------*/
        /** \brief Constructor of the edge going from \p ATail to \p AHead
         *
         *  \param[in] ATail the edge tail
         *  \param[in] AHead the edge head
         */
        /*--------------------------------------------------------------------*/
        GraphEdge(const TCellID AID, GraphNode* ATail=0,  GraphNode* AHead=0);
        /*--------------------------------------------------------------------*/
        /** \brief Returns the edge id
         */
        /*--------------------------------------------------------------------*/
        TCellID id();

        /*--------------------------------------------------------------------*/
        /** \brief Access to the tail node
         *
         *  \return an pointer on the tail node
         */
        /*--------------------------------------------------------------------*/
        GraphNode* tail();
        /*--------------------------------------------------------------------*/
        /** \brief Access to the head node
         *
         *  \return an pointer on the head node
         */
        /*--------------------------------------------------------------------*/
        GraphNode* head();
        
        std::vector<GraphEdge*> getEdgesStartingFrom(GraphNode* ANode);
    private:
        TCellID m_id;
        GraphNode* m_tail;
        GraphNode* m_head;
    };
    
    /*------------------------------------------------------------------------*/
    /** \class OrientedGraph
     *
     *  \brief Class defining a simple oriented graph. Useful for some pure
     *         topoligical algorithm on meshes.
     */
    /*------------------------------------------------------------------------*/
    class OrientedGraph{
        
    public:
        /*--------------------------------------------------------------------*/
        /** \brief Constructor
         *
         * \param ANbNodes the number of nodes in the graph. Nodes will be 
         *        numbered from 0 to ANbNodes-1
         */
        /*--------------------------------------------------------------------*/
        OrientedGraph(const int ANbNodes);
        /*--------------------------------------------------------------------*/
        /** \brief Constructor
         *
         * \param AIDS the ids we want to put on each nodes
         */
        /*--------------------------------------------------------------------*/
        OrientedGraph(const std::set<TCellID>& AIDS);
        /*--------------------------------------------------------------------*/
        /** \brief Destructor
         */
        /*--------------------------------------------------------------------*/
        virtual ~OrientedGraph();
        
        /*--------------------------------------------------------------------*/
        /** \brief Add the oriented edge [\p ATail, \p AHead]. Warning, be 
         *         careful to use graph node stored in this graph. No check 
         *         done.
         *
         *  \param[in] AID  the edge id
         *  \param[in] ATail the edge tail
         *  \param[in] AHead the edge head
         */
        /*--------------------------------------------------------------------*/
        bool addEdge(const TCellID AID,
                     const TCellID ATail,
                     const TCellID AHead);
        
        /*--------------------------------------------------------------------*/
        /** \brief Returs a reference on the node indexed \p AIndex
         *
         *  \return a node reference
         */
        /*--------------------------------------------------------------------*/
        GraphNode* node(const TCellID AIndex);
        /*--------------------------------------------------------------------*/
        /** \brief Returs a reference on the edge indexed \p AIndex
         *
         *  \return an edge reference
         */
        /*--------------------------------------------------------------------*/
        GraphEdge* edge(const TCellID AIndex);
        /*--------------------------------------------------------------------*/
        /** \brief Returns the number of nodes
         *
         *  \return the number of nodes
         */
        /*--------------------------------------------------------------------*/
        int nbNodes() const;
        /*--------------------------------------------------------------------*/
        /** \brief Returns the number of edges
         *
         *  \return the number of edges
         */
        /*--------------------------------------------------------------------*/
        int nbEdges() const;
        
        void updateNodes();
        friend std::ostream& operator<<(std::ostream& AStr, const OrientedGraph& AG);
    
        
        std::vector<GraphNode> m_nodes;
        std::map<TCellID, GraphNode*> m_map_nodes;
        std::vector<GraphEdge> m_edges;
    };
    
    
    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_ORIENTED_GRAPH_H_ */
/*----------------------------------------------------------------------------*/
