#include "gmds/medialaxis/Conformalizer.h"
#include <queue>
#include <time.h> 
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
Conformalizer::Conformalizer(gmds::Mesh &AMesh, std::vector<NonConformalHalfEdge> AHalfEdges, std::vector<int> ALengths)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
    // Old node/new Node correspondance
    m_mesh->newVariable<int,GMDS_NODE>("node2group");
    // List of ordered new nodes of each old edge
    m_mesh->newVariable<std::vector<Node>,GMDS_EDGE>("oldEdges2newNodes");
    // Old faces/new nodes correspondance
    m_mesh->newVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2newNodes");
    // Internal constraints
    try {m_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");}
    catch (GMDSException& e){m_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");}
    try {m_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");}
    catch (GMDSException& e){m_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");}
    // Correspondance with min tri nodes
	try {m_mesh->getVariable<int,GMDS_NODE>("TMeshNode2MinTriNode");}
    catch (GMDSException& e){m_mesh->newVariable<int,GMDS_NODE>("TMeshNode2MinTriNode");}
    // Corners
	try {m_mesh->getVariable<int,GMDS_NODE>("corner");}
    catch (GMDSException& e){m_mesh->newVariable<int,GMDS_NODE>("corner");}

	// Corresponding non-conformal half edges
	m_half_edges = AHalfEdges;
	// Half edges length. WARNING : all lengths must be > 0
	m_half_edges_lengths = ALengths;

	// Conformal mesh
	m_conformal_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
    // Node/group of intermediate nodes correspondance
    m_conformal_mesh->newVariable<int,GMDS_NODE>("node2group");
    // Node/representative of the group correspondance
    m_conformal_mesh->newVariable<int,GMDS_NODE>("node2representative");
    // Mark with 1 edges corresponding to internal constraints
	m_conformal_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");
	// Mark with 1 nodes belonging to an internal constraint
	m_conformal_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    // Boundary & constraint nodes/reference node correspondance
    m_conformal_mesh->newVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    // Mark with 1 corners of the boundary and internal constraints
    m_conformal_mesh->newVariable<int,GMDS_NODE>("corner");

    // Intermediate mesh                                          
	m_intermediate_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
    m_intermediate_mesh->newVariable<int,GMDS_EDGE>("length");
}

/*----------------------------------------------------------------------------*/
void Conformalizer::writeConformalMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the conformal mesh"<<std::endl;
	IGMeshIOService ioService(m_conformal_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void Conformalizer::writeIntermediateMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the intermediate mesh"<<std::endl;
	IGMeshIOService ioService(m_intermediate_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
Mesh Conformalizer::getConformalMesh()
{
	return *m_conformal_mesh;
}

/*----------------------------------------------------------------------------*/
void Conformalizer::setConformalMeshConnectivity()
{
	std::cout<<"> Setting conformal mesh connectivity"<<std::endl;
	MeshDoctor doc(m_conformal_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void Conformalizer::setIntermediateMeshConnectivity()
{
	std::cout<<"> Setting intermediate mesh connectivity"<<std::endl;
	MeshDoctor doc(m_intermediate_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildNonConformalNodesGroups()
{
    std::cout<<"> Building groups of nodes of the non conformal mesh"<<std::endl;
    auto n2g = m_mesh->getVariable<int,GMDS_NODE>("node2group");
    auto l = m_mesh->getVariable<int,GMDS_EDGE>("length");
    int groupId = 0;
    std::vector<std::vector<TCellID>> groups;
    std::vector<TCellID> newGroup;
    auto visited = m_mesh->newMark<Node>();
    for (auto n_id:m_mesh->nodes())
    {
        if (!m_mesh->isMarked<Node>(n_id,visited))
        {
            // We create a new group
            newGroup.clear();
            std::queue<Node> toAdd;
            Node n = m_mesh->get<Node>(n_id);
            toAdd.push(n);
            m_mesh->mark<Node>(n_id,visited);
            while (!toAdd.empty())
            {
                n = toAdd.front();
                newGroup.push_back(n.id());
                n2g->set(n.id(),groupId);
                toAdd.pop();
                for (auto e:n.get<Edge>())
                {
                    if (l->value(e.id()) == 0)
                    {
                        for (auto n1:e.get<Node>())
                        {
                            if (!m_mesh->isMarked<Node>(n1.id(),visited))
                            {
                                toAdd.push(n1);
                                m_mesh->mark<Node>(n1.id(),visited);
                            }
                        }
                    }
                }
            }
            groups.push_back(newGroup);
            groupId += 1;
        }
    }
    m_non_conformal_nodes_groups = groups;
}

/*----------------------------------------------------------------------------*/
void Conformalizer::addOldNodesOnConformalMesh()
{
    std::cout<<"> Adding old nodes on the conformal mesh"<<std::endl;
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    auto tmn2mtn = m_mesh->getVariable<int,GMDS_NODE>("TMeshNode2MinTriNode");
    auto constr = m_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto c1 = m_mesh->getVariable<int,GMDS_NODE>("corner");
    auto c2 = m_conformal_mesh->getVariable<int,GMDS_NODE>("corner");
    for (auto group:m_non_conformal_nodes_groups)
    {
        double x = 0.;
        double y = 0.;
        Node rep;
        bool boundOrConstr = false;
        bool cor = false;
        for (auto n_id:group)
        {
            Node n = m_mesh->get<Node>(n_id);
            x += n.X();
            y += n.Y();
            if (!isInterior(n) || constr->value(n_id) == 1)
            {
                boundOrConstr = true;
                rep = n;
            }
            if (c1->value(n_id) == 1)
                cor = true;
        }
        x = x/double(group.size());
        y = y/double(group.size());
        math::Point P(x,y,0.);
        Node new_node = m_conformal_mesh->newNode(P);
        if (boundOrConstr)
            n2rn->set(new_node.id(),tmn2mtn->value(rep.id()));
        else
            n2rn->set(new_node.id(),-1);
        if (cor)
            c2->set(new_node.id(),1);
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::addSubdividingNodes()
{
    std::cout<<"> Adding subdividing nodes on the conformal mesh"<<std::endl;
    auto l = m_mesh->getVariable<int,GMDS_EDGE>("length");
    auto oe2nn = m_mesh->getVariable<std::vector<Node>,GMDS_EDGE>("oldEdges2newNodes");
    auto n2g = m_mesh->getVariable<int,GMDS_NODE>("node2group");
    auto int_constr = m_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto constr_nodes = m_conformal_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    // Add nodes subdividing edges
    for (auto e_id:m_mesh->edges())
    {
        int len = l->value(e_id);
        Edge e = m_mesh->get<Edge>(e_id);
        Node n1 = e.get<Node>()[0];
        Node n2 = e.get<Node>()[1];
        bool boundOrConstr = (n2rn->value(n1.id()) >= 0 && n2rn->value(n2.id()) >= 0);
        std::vector<Node> ordered_nodes;
        if (len == 0)
        {
            ordered_nodes.push_back(m_conformal_mesh->get<Node>(n2g->value(n1.id())));
            oe2nn->set(e_id,ordered_nodes);
        }
        if (len >= 1)
        {
            ordered_nodes.push_back(m_conformal_mesh->get<Node>(n2g->value(n1.id())));
            math::Point dir = n2.point()+(-1.)*n1.point();
            for (int i = 1; i < len; i++)
            {
                Node newNode = m_conformal_mesh->newNode(n1.point()+(double(i)/double(len))*dir);
                ordered_nodes.push_back(newNode);
                if (boundOrConstr)
                    n2rn->set(newNode.id(),n2rn->value(n1.id()));
                else
                    n2rn->set(newNode.id(),-1);
            }
            ordered_nodes.push_back(m_conformal_mesh->get<Node>(n2g->value(n2.id())));
            oe2nn->set(e_id,ordered_nodes);
        }
    }
    // Mark the nodes of the conformal mesh belonging to constraints
    for (auto e_id:m_mesh->edges())
    {
        if (int_constr->value(e_id) == 1)
        {
            for (auto n:oe2nn->value(e_id))
            {
                constr_nodes->set(n.id(),1);
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildHalfEdges2newNodesConnectivity()
{
    std::cout<<"> Building 'half edge 2 new nodes' connectivity"<<std::endl;
    auto oe2nn = m_mesh->getVariable<std::vector<Node>,GMDS_EDGE>("oldEdges2newNodes");
    auto n2g = m_mesh->getVariable<int,GMDS_NODE>("node2group");
    std::vector<std::vector<TCellID>> he2nn(m_half_edges.size());
    for (auto he:m_half_edges)
    {
        std::vector<Node> ordered_old_nodes = he.getOrderedNodes();
        std::vector<TCellID> ordered_new_nodes;
        ordered_new_nodes.push_back(n2g->value(he.firstNode().id()));
        int previous_id = n2g->value(he.firstNode().id());
        for (int i = 0; i < ordered_old_nodes.size()-1; i++)
        {
            Node n1 = ordered_old_nodes[i];
            Node n2 = ordered_old_nodes[i+1];
            Edge e = getEdge(n1,n2);
            std::vector<Node> nodes = oe2nn->value(e.id());
            if (e.get<Node>()[1].id() == n1.id())
                std::reverse(nodes.begin(),nodes.end());
            for (int j = 1; j < nodes.size(); j++)
            {
                if (previous_id != nodes[j].id())
                {
                    ordered_new_nodes.push_back(nodes[j].id());
                    previous_id = nodes[j].id();
                } 
            }
        }
        he2nn[he.id()] = ordered_new_nodes;
    }
    m_half_edges_to_new_nodes = he2nn;
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildIntermediateMeshNodes()
{
    std::cout<<"> Building intermediate mesh nodes"<<std::endl;
    for (auto n_id:m_conformal_mesh->nodes())
        m_intermediate_mesh->newNode(m_conformal_mesh->get<Node>(n_id).point());
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildIntermediateMeshEdges()
{
    std::cout<<"> Building intermediate mesh edges"<<std::endl;
    auto l = m_intermediate_mesh->getVariable<int,GMDS_EDGE>("length");
    // Build edges between subdivision nodes
    for (auto f_id:m_mesh->faces())
    {
        // Get the half edges of the face
		NonConformalHalfEdge e1 = m_half_edges[4*f_id];
		NonConformalHalfEdge e2 = m_half_edges[4*f_id+1];
		NonConformalHalfEdge e3 = m_half_edges[4*f_id+2];
		NonConformalHalfEdge e4 = m_half_edges[4*f_id+3];
        std::vector<TCellID> nodes1 = m_half_edges_to_new_nodes[e1.id()];
        std::vector<TCellID> nodes2 = m_half_edges_to_new_nodes[e2.id()];
        std::vector<TCellID> nodes3 = m_half_edges_to_new_nodes[e3.id()];
        std::vector<TCellID> nodes4 = m_half_edges_to_new_nodes[e4.id()];
        int I = nodes1.size();
        int J = nodes2.size();
        for (int i = 1; i < I-1; i++)
            {
                Edge e = m_intermediate_mesh->newEdge(nodes1[i],nodes3[I-1-i]);
                l->set(e.id(),m_half_edges_lengths[e2.id()]);
            }
        for (int i = 1; i < J-1; i++)
            {
                Edge e = m_intermediate_mesh->newEdge(nodes2[i],nodes4[J-1-i]);
                l->set(e.id(),m_half_edges_lengths[e3.id()]);
            }
    }
    
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildIntermediateNodesGroups()
{
    std::cout<<"> Building intermediate nodes groups"<<std::endl;
    auto n2g = m_conformal_mesh->getVariable<int,GMDS_NODE>("node2group");
    auto l = m_intermediate_mesh->getVariable<int,GMDS_EDGE>("length");
    int groupId = 0;
    std::vector<std::vector<TCellID>> groups;
    std::vector<TCellID> newGroup;
    auto visited = m_intermediate_mesh->newMark<Node>();
    for (auto n_id:m_intermediate_mesh->nodes())
    {
        if (!m_intermediate_mesh->isMarked<Node>(n_id,visited))
        {
            // We create a new group
            newGroup.clear();
            std::queue<Node> toAdd;
            Node n = m_intermediate_mesh->get<Node>(n_id);
            toAdd.push(n);
            m_intermediate_mesh->mark<Node>(n_id,visited);
            while (!toAdd.empty())
            {
                n = toAdd.front();
                newGroup.push_back(n.id());
                n2g->set(n.id(),groupId);
                toAdd.pop();
                for (auto e:n.get<Edge>())
                {
                    if (l->value(e.id()) == 0)
                    {
                        for (auto n1:e.get<Node>())
                        {
                            if (!m_intermediate_mesh->isMarked<Node>(n1.id(),visited))
                            {
                                toAdd.push(n1);
                                m_intermediate_mesh->mark<Node>(n1.id(),visited);
                            }
                        }
                    }
                }
            }
            groups.push_back(newGroup);
            groupId += 1;
        }
    }
    m_intermediate_nodes_groups = groups;
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildRepresentativeNodes()
{
    std::cout<<"> Building representative nodes"<<std::endl;
    auto constr = m_conformal_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto n2r = m_conformal_mesh->getVariable<int,GMDS_NODE>("node2representative");
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    for (auto group:m_intermediate_nodes_groups)
    {
        if (group.size() == 1)
            n2r->set(group[0],group[0]);
        else
        {
            double x = 0.;
            double y = 0.;
            bool constrained = false;
            for (auto n_id:group)
            {
                x += m_conformal_mesh->get<Node>(n_id).X();
                y += m_conformal_mesh->get<Node>(n_id).Y();
                if (constr->value(n_id) == 1)
                    constrained = true;
            }
            x = x/double(group.size());
            y = y/double(group.size());
            math::Point P(x,y,0.);
            Node n = m_conformal_mesh->newNode(P);
            n2r->set(n.id(),n.id());
            for (auto n_id:group)
                n2r->set(n_id,n.id());
            if (constrained)
                constr->set(n.id(),1);
            n2rn->set(n.id(),-1);
        }
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::builIntNodesAndFaces2newNodesConnectivity()
{
    std::cout<<"> Building nodes of the conformal mesh internal to faces of the non conformal mesh"<<std::endl;
    auto f2nn = m_mesh->getVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2newNodes");
    auto n2r = m_conformal_mesh->getVariable<int,GMDS_NODE>("node2representative");
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
	for (auto f_id:m_mesh->faces())
	{
		// Get the half edges of the face
		NonConformalHalfEdge e1 = m_half_edges[4*f_id];
		NonConformalHalfEdge e2 = m_half_edges[4*f_id+1];
		NonConformalHalfEdge e3 = m_half_edges[4*f_id+2];
		NonConformalHalfEdge e4 = m_half_edges[4*f_id+3];
		std::vector<std::vector<int>> nodesIDs(m_half_edges_lengths[e2.id()]+1);
		std::vector<int> row(m_half_edges_lengths[e1.id()]+1);
		// First row
		for (int i = 0; i < m_half_edges_lengths[e1.id()]+1; i++)
			row[i] = n2r->value(m_half_edges_to_new_nodes[e1.id()][i]);
		nodesIDs[0] = row;
		// Last row
		for (int i = 0; i < m_half_edges_lengths[e3.id()]+1; i++)
			row[i] = n2r->value(m_half_edges_to_new_nodes[e3.id()][m_half_edges_lengths[e3.id()]-i]);
		nodesIDs[m_half_edges_lengths[e2.id()]] = row;
		for (int j = 1; j < m_half_edges_lengths[e2.id()]; j++)
		{
			row[0] = n2r->value(m_half_edges_to_new_nodes[e4.id()][m_half_edges_lengths[e4.id()]-j]);
			row[m_half_edges_lengths[e1.id()]] = n2r->value(m_half_edges_to_new_nodes[e2.id()][j]);
			for (int i = 1; i < m_half_edges_lengths[e1.id()]; i++)
			{
				math::Point P1 = m_conformal_mesh->get<Node>(n2r->value(m_half_edges_to_new_nodes[e1.id()][i])).point();
				math::Point P2 = m_conformal_mesh->get<Node>(n2r->value(m_half_edges_to_new_nodes[e3.id()][m_half_edges_lengths[e3.id()]-i])).point();
				math::Point P = P1 + (double(j)/double(m_half_edges_lengths[e2.id()]))*(P2+(-1.)*P1);
				Node n = m_conformal_mesh->newNode(P);
				row[i] = n.id();
                n2r->set(n.id(),n.id());
                n2rn->set(n.id(),-1);
			}
			nodesIDs[j] = row;
		}
		f2nn->set(f_id,nodesIDs);
	}
}

/*----------------------------------------------------------------------------*/
void Conformalizer::buildConformalMeshFaces()
{
	std::cout<<"> Building faces of the conformal mesh"<<std::endl;
	auto f2nn = m_mesh->getVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2newNodes");
	for (auto f_id:m_mesh->faces())
	{
		std::vector<std::vector<int>> nodes = f2nn->value(f_id);
		int I = nodes.size();
		int J = nodes[0].size();
		TCellID i1,i2,i3,i4;
		for (int i = 0; i < I-1; i++)
		{
			for (int j = 0; j < J-1; j++)
			{
				i1 = nodes[i][j];
				i2 = nodes[i+1][j];
				i3 = nodes[i+1][j+1];
				i4 = nodes[i][j+1];
				m_conformal_mesh->newFace({i1,i2,i3,i4});
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void Conformalizer::deleteSuperfluousNodes()
{
    std::cout<<"> Deleting superfluous nodes of the conformal mesh"<<std::endl;
    auto n2r = m_conformal_mesh->getVariable<int,GMDS_NODE>("node2representative");
    for (auto n_id:m_conformal_mesh->nodes())
    {
        if (n2r->value(n_id) != n_id)
            m_conformal_mesh->deleteNode(n_id);
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::markInternalConstraints()
{
    std::cout<<"> Marking internal constraints on the conformal mesh"<<std::endl;
    auto int_constr = m_conformal_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
	auto nodes_constr = m_conformal_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    for (auto e_id:m_conformal_mesh->edges())
    {
		Edge e = m_conformal_mesh->get<Edge>(e_id);
		Node n1 = e.get<Node>()[0];
		Node n2 = e.get<Node>()[1];
		if (nodes_constr->value(n1.id()) == 1 && nodes_constr->value(n2.id()) == 1)
			int_constr->set(e.id(),1);
		else if (nodes_constr->value(n1.id()) == 1 && isInterior(n1) && !isInterior(n2))
			int_constr->set(e.id(),1);	
		else if (nodes_constr->value(n2.id()) == 1 && isInterior(n2) && !isInterior(n1))
			int_constr->set(e.id(),1);	
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::execute()
{
    std::cout<<" "<<std::endl;
    std::cout<<"=== Building a conformal mesh from from a non conformal mesh and quantization results ==="<<std::endl;
    buildNonConformalNodesGroups();
    addOldNodesOnConformalMesh();
    addSubdividingNodes();
    buildHalfEdges2newNodesConnectivity();
    buildIntermediateMeshNodes();
    buildIntermediateMeshEdges();
    setIntermediateMeshConnectivity();
    buildIntermediateNodesGroups();
    buildRepresentativeNodes();
    builIntNodesAndFaces2newNodesConnectivity();
    buildConformalMeshFaces();
    deleteSuperfluousNodes();
    setConformalMeshConnectivity();
    markInternalConstraints();
    std::cout<<"========================================================================================"<<std::endl;
    std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
void Conformalizer::projectOnBoundary(Mesh &ARefMesh)
{
    std::cout<<"> Projecting on the boundary and the internal constraints"<<std::endl;
    auto nodes_constr = m_conformal_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    for (auto n_id:m_conformal_mesh->nodes())
    {
        Node n = m_conformal_mesh->get<Node>(n_id);
        if (!isInterior(n) || nodes_constr->value(n_id) == 1)
        {
            double min_dist = 1e6;
            Node closest;
            for (auto boundNode_id:ARefMesh.nodes())
            {
                Node boundNode = ARefMesh.get<Node>(boundNode_id);
                if (boundNode.point().distance(n.point()) < min_dist)
                {
                    min_dist = boundNode.point().distance(n.point());
                    closest = boundNode;
                }
            }
            n.setX(closest.X());
            n.setY(closest.Y());
            n2rn->set(n.id(),closest.id());
        }
    }
}

/*----------------------------------------------------------------------------*/
void Conformalizer::smooth(Mesh &ARefMesh)
{
    std::cout<<"> Smoothing the conformal mesh"<<std::endl;
    auto nodes_constr = m_conformal_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto int_constr = m_conformal_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto n2rn = m_conformal_mesh->getVariable<int,GMDS_NODE>("Boundary&ConstraintNode2ReferenceNode"); 
    auto corner = m_conformal_mesh->getVariable<int,GMDS_NODE>("corner");
    for (auto n_id:m_conformal_mesh->nodes())
    {
        Node n = m_conformal_mesh->get<Node>(n_id);
        if (isInterior(n))
        {
            if (nodes_constr->value(n.id()) == 0)
            {
                double x = 0.;
                double y = 0.;
                double NbNeighbours = 0.;
                for (auto e:n.get<Edge>())
                {
                    for (auto n1:e.get<Node>())
                    {
                        if (n1.id() != n.id())
                        {
                            x += n1.X();
                            y += n1.Y();
                            NbNeighbours += 1.;
                        }
                    }
                }
                x = x/NbNeighbours;
                y = y/NbNeighbours;
                n.setX(x);
                n.setY(y);
            }
            else
            {
                int NbConstraintedEdges = 0;
                for (auto e:n.get<Edge>())
                {
                    if (int_constr->value(e.id()) == 1)
                        NbConstraintedEdges += 1;
                }
                if (NbConstraintedEdges == 2 && corner->value(n_id) == 0)
                {
                    // Find the two constrainted edges
                    Edge e1,e2;
                    for (auto e:n.get<Edge>())
                    {
                        if (int_constr->value(e.id()) == 1)
                        {
                            e1 = e;
                            break;
                        }
                    }
                    for (auto e:n.get<Edge>())
                    {
                        if (int_constr->value(e.id()) == 1 && e.id() != e1.id())
                        {
                            e2 = e;
                            break;
                        }
                    }


                    double x = 0.;
                    double y = 0.;
                    
                    double NbNeighbours = 0.;
                    for (auto e:n.get<Edge>())
                    {
                        for (auto n1:e.get<Node>())
                        {
                            if (n1.id() != n.id())
                            {
                                x += n1.X();
                                y += n1.Y();
                                NbNeighbours += 1.;
                            }
                        }
                    }
                    x = x/NbNeighbours;
                    y = y/NbNeighbours;
                    math::Point P(x,y,0.);

                    // Projection on the boundary using the reference mesh
                    Node n1,n2;
                    if (e1.get<Node>()[0].id() == n.id())
                        n1 = e1.get<Node>()[1];
                    else
                        n1 = e1.get<Node>()[0];
                    if (e2.get<Node>()[0].id() == n.id())
                        n2 = e2.get<Node>()[1];
                    else
                        n2 = e2.get<Node>()[0];

                    Node boundNode1 = ARefMesh.get<Node>(n2rn->value(n1.id()));
                    Node boundNode2 = ARefMesh.get<Node>(n2rn->value(n2.id()));
                    std::vector<Node> arc = shortestPathAlongBoundaryOrConstraints(boundNode1,boundNode2,ARefMesh);
                    double min_dist = 1e6;
                    Node closest;
                    for (auto boundNode:arc)
                    {
                        if (boundNode.point().distance(P) < min_dist)
                        {
                            min_dist = boundNode.point().distance(P);
                            closest = boundNode;
                        }
                    }
                    // for (auto boundNode_id:ARefMesh.nodes())
                    // {
                    //     Node boundNode = ARefMesh.get<Node>(boundNode_id);
                    //     if (boundNode.point().distance(P) < min_dist)
                    //     {
                    //         min_dist = boundNode.point().distance(P);
                    //         closest = boundNode;
                    //     }
                    // }
                    n.setX(closest.X());
                    n.setY(closest.Y());
                    n2rn->set(n.id(),closest.id());
                }
            }
        }
        else
        {
            if (nodes_constr->value(n.id()) == 0 && corner->value(n.id()) == 0) // Prevent the corners from moving
            //if (nodes_constr->value(n.id()) == 0 ) // Authorize the corners to move
            {
                // Find the two boundary edges containing n
                Edge e1,e2;
                for (auto e:n.get<Edge>())
                {
                    if (e.get<Face>().size() == 1)
                    {
                        e1 = e;
                        break;
                    }
                }
                for (auto e:n.get<Edge>())
                {
                    if (e.get<Face>().size() == 1 && e.id() != e1.id())
                    {
                        e2 = e;
                        break;
                    }
                }

                double x = 0.;
                double y = 0.;
                
                double NbNeighbours = 0.;
                for (auto e:n.get<Edge>())
                {
                    for (auto n1:e.get<Node>())
                    {
                        if (n1.id() != n.id())
                        {
                            x += n1.X();
                            y += n1.Y();
                            NbNeighbours += 1.;
                        }
                    }
                }
                x = x/NbNeighbours;
                y = y/NbNeighbours;
                math::Point P(x,y,0.);

                // Projection on the boundary using the reference mesh
                Node n1,n2;
                if (e1.get<Node>()[0].id() == n.id())
                    n1 = e1.get<Node>()[1];
                else
                    n1 = e1.get<Node>()[0];
                if (e2.get<Node>()[0].id() == n.id())
                    n2 = e2.get<Node>()[1];
                else
                    n2 = e2.get<Node>()[0];

                Node boundNode1 = ARefMesh.get<Node>(n2rn->value(n1.id()));
                Node boundNode2 = ARefMesh.get<Node>(n2rn->value(n2.id()));
                std::vector<Node> arc = shortestPathAlongBoundaryOrConstraints(boundNode1,boundNode2,ARefMesh);
                double min_dist = 1e6;
                Node closest;
                for (auto boundNode:arc)
                {
                    if (boundNode.point().distance(P) < min_dist)
                    {
                        min_dist = boundNode.point().distance(P);
                        closest = boundNode;
                    }
                }
                // for (auto boundNode_id:ARefMesh.nodes())
                // {
                //     Node boundNode = ARefMesh.get<Node>(boundNode_id);
                //     if (boundNode.point().distance(P) < min_dist)
                //     {
                //         min_dist = boundNode.point().distance(P);
                //         closest = boundNode;
                //     }
                // }
                n.setX(closest.X());
                n.setY(closest.Y());
                n2rn->set(n.id(),closest.id());
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/