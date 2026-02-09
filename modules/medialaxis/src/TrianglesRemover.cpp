#include "gmds/medialaxis/TrianglesRemover.h"
#include <queue>
#include <time.h> 
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
TrianglesRemover::TrianglesRemover(gmds::Mesh &AMesh)
{
	// Non conformal quad mesh
	m_degenerated_mesh = &AMesh;
    // Degenerated mesh nodes/T-mesh nodes correspondance
    m_degenerated_mesh->newVariable<int,GMDS_NODE>("degMeshNode2TMeshNode");
    // Node/group correspondance
    m_degenerated_mesh->newVariable<int,GMDS_NODE>("degenerated_node_to_group");
    
    // Triangle fan id
    try {m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_NODE>("fan_id");}
    // Triangle fans angle
    try {m_degenerated_mesh->getVariable<double,GMDS_NODE>("fan_angle");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<double,GMDS_NODE>("fan_angle");}
    // Separatrices
    try {m_degenerated_mesh->getVariable<int,GMDS_EDGE>("belongs_to_a_separatrix");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_EDGE>("belongs_to_a_separatrix");}
    // Corners
	try {m_degenerated_mesh->getVariable<int,GMDS_NODE>("corner");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_NODE>("corner");}
    // Internal constraints
    try {m_degenerated_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");}
    try {m_degenerated_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");}
    // Triangles
    try {m_degenerated_mesh->getVariable<int,GMDS_FACE>("is_a_triangle");}
    catch (GMDSException& e){m_degenerated_mesh->newVariable<int,GMDS_FACE>("is_a_triangle");}

    // T-mesh                                          
	m_t_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
    // T-junctions on the T-mesh
    m_t_mesh->newVariable<int,GMDS_NODE>("is_a_T-junction");
    // Corners
    m_t_mesh->newVariable<int,GMDS_NODE>("corner");
    // Internal constraints
    m_t_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");
    m_t_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    // Marking with 1 triangles created by remeshing the triangles fans
    m_t_mesh->newVariable<int,GMDS_FACE>("is_an_unwanted_triangle");
    // Marking with 1 vertices of unwanted triangles
    m_t_mesh->newVariable<int,GMDS_NODE>("is_an_unwanted_triangle_vertex");
    // Marking with 1 nodes created by removing the unwanted triangles
    m_t_mesh->newVariable<int,GMDS_NODE>("is_a_new_node");

    // Final T-mesh                                          
	m_final_t_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
    // T-junctions on the T-mesh
    m_final_t_mesh->newVariable<int,GMDS_NODE>("is_a_T-junction");
    // Corners
    m_final_t_mesh->newVariable<int,GMDS_NODE>("corner");
    // Internal constraints
    m_final_t_mesh->newVariable<int,GMDS_EDGE>("internal_constraint");
    m_final_t_mesh->newVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    // Marking with 1 nodes created by removing the unwanted triangles
    m_final_t_mesh->newVariable<int,GMDS_NODE>("is_a_new_node");
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::writeDegeneratedMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the degenerated mesh"<<std::endl;
	IGMeshIOService ioService(m_degenerated_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::writeTMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the T-mesh"<<std::endl;
	IGMeshIOService ioService(m_t_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
std::vector<Face> TrianglesRemover::fan(int AFanId)
{
    auto fan_id = m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");
    auto separatrix = m_degenerated_mesh->getVariable<int,GMDS_EDGE>("belongs_to_a_separatrix");
    auto tri = m_degenerated_mesh->getVariable<int,GMDS_FACE>("is_a_triangle");
    std::vector<Face> output_fan;
    std::queue<Face> faces_to_add;
    auto open = m_degenerated_mesh->newVariable<int,GMDS_EDGE>("open");
    auto visited = m_degenerated_mesh->newVariable<int,GMDS_FACE>("visited");
    // Find the origin triangles of the fan

    for (auto f_id:m_degenerated_mesh->faces())
    {
        if (tri->value(f_id) == 1)
        {
            Face f = m_degenerated_mesh->get<Face>(f_id);
            int nb_vertices = 0;
            Node n1,n2;
            bool n1_found = false;
            for (auto n:f.get<Node>())
            {
                if (fan_id->value(n.id()) == AFanId)
                {
                    nb_vertices += 1;
                    if (n1_found)
                        n2 = n;
                    else
                    {
                        n1 = n;
                        n1_found = true;
                    }
                }  
            }
            if (nb_vertices == 2)
            {
                faces_to_add.push(f);
                visited->set(f_id,1);
                Edge e1 = getEdge(n1,n2);
                Edge e2 = opp(f,e1);
                if (separatrix->value(e2.id()) == 0)
                    open->set(e2.id(),1);
            }
        }
    }

    // Fill the fan
    while (!faces_to_add.empty())
    {
        Face f1 = faces_to_add.front();
        faces_to_add.pop();
        output_fan.push_back(f1);
        for (auto e:f1.get<Edge>())
        {
            if (open->value(e.id()) == 1)
            {
                for (auto f2:e.get<Face>())
                {
                    if (visited->value(f2.id()) == 0)
                    {
                        faces_to_add.push(f2);
                        visited->set(f2.id(),1);
                        Edge e2 = opp(f2,e);
                        if (separatrix->value(e2.id()) == 0)
                            open->set(e2.id(),1);
                    }
                }
            }
        }
    }

    // // To visualize the fan
    // auto bttf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("belongs_to_the_fan");
    // for (auto f:output_fan)
    //     bttf->set(f.id(),1);

    m_degenerated_mesh->deleteVariable(GMDS_EDGE,open);
    m_degenerated_mesh->deleteVariable(GMDS_FACE,visited);
    return output_fan;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> TrianglesRemover::outline(std::vector<Face> &AV)
{
    auto n2g = m_degenerated_mesh->getVariable<int,GMDS_NODE>("degenerated_node_to_group");
    // Mark the faces of the shape
    auto bttf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("belongs_to_the_fan");
    for (auto f:AV)
        bttf->set(f.id(),1);
    // Mark the edges of the outline
    auto outline_edges = m_degenerated_mesh->newVariable<int,GMDS_EDGE>("outline");
    for (auto f:AV)
    {
        for (auto e:f.get<Edge>())
        {
            if (e.get<Face>().size() == 1)
                outline_edges->set(e.id(),1);
            for (auto f2:e.get<Face>())
            {
                if (bttf->value(f2.id()) == 0)
                    outline_edges->set(e.id(),1);
            }
        }
    }
    // Find the nodes of the outline
    std::vector<Node> outline_nodes;
    // Find a face and an edge at the border of the fan
    Face f_start;
    Edge e_start;
    bool found = false;
    for (auto f:AV)
    {
        for (auto e:f.get<Edge>())
        {
            if (outline_edges->value(e.id()) == 1)
            {
                if (n2g->value(e.get<Node>()[0].id()) != n2g->value(e.get<Node>()[1].id()))
                {
                    f_start = f;
                    e_start = e;
                    found = true;
                    break;
                }
                
            }
        }
        if (found)
            break;
    }
    // // Orientate e_start with respect to f_start
    // int pos1,pos2;
    // for (int i = 0; i < 4; i++)
    // {
    //     if (e_start.get<Node>()[0].id() == f_start.get<Node>()[i].id())
    //         pos1 = i;
    //     if (e_start.get<Node>()[1].id() == f_start.get<Node>()[i].id())
    //         pos2 = i;
    // }
    // Node prev,nxt;
    // if ((pos1 < 3 && pos1 < pos2) || (pos1 == 3 && pos2 == 0))
    // {
    //     prev = e_start.get<Node>()[0];
    //     nxt = e_start.get<Node>()[1];
    // }
    // else
    // {
    //     prev = e_start.get<Node>()[1];
    //     nxt = e_start.get<Node>()[0];
    // }
    Node prev = m_degenerated_nodes_groups[n2g->value(e_start.get<Node>()[0].id())][0];
    Node nxt = m_degenerated_nodes_groups[n2g->value(e_start.get<Node>()[1].id())][0];

    outline_nodes.push_back(prev);
    Node old_prev, candidate;
    while (nxt.id() != outline_nodes[0].id())
    {
        old_prev = prev;
        prev = nxt;
        outline_nodes.push_back(prev);
        bool found = false;
        for (auto n:m_degenerated_nodes_groups[n2g->value(prev.id())])
        {
            for (auto e:n.get<Edge>())
            {
                if (n2g->value(e.get<Node>()[0].id()) == n2g->value(e.get<Node>()[1].id()))
                    continue;
                if (n2g->value(e.get<Node>()[0].id()) == n2g->value(prev.id()))
                    candidate = m_degenerated_nodes_groups[n2g->value(e.get<Node>()[1].id())][0];
                else
                    candidate = m_degenerated_nodes_groups[n2g->value(e.get<Node>()[0].id())][0];
                if (outline_edges->value(e.id()) == 1 && candidate.id() != old_prev.id() && n2g->value(candidate.id()) != n2g->value(n.id()))
                {
                    nxt = candidate;
                    found = true;
                    break;
                }
            }
            if (found)
                break;
        }
        // for (auto e:prev.get<Edge>())
        // {
        //     if (e.get<Node>()[0].id() == prev.id())
        //         candidate = e.get<Node>()[1];
        //     else
        //         candidate = e.get<Node>()[0];
        //     if (outline_edges->value(e.id()) == 1 && candidate.id() != old_prev.id())
        //     {
        //         nxt = candidate;
        //         break;
        //     }
        // }
    }


    // // To visualize the outline
    // auto btto = m_degenerated_mesh->newVariable<int,GMDS_NODE>("belongs_to_the_outline");
    // for (auto n:outline_nodes)
    //     btto->set(n.id(),1);

    m_degenerated_mesh->deleteVariable(GMDS_FACE,bttf);
    m_degenerated_mesh->deleteVariable(GMDS_EDGE,outline_edges);

    return outline_nodes;
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::buildTMesh()
{
    std::cout<<"> Building the T-mesh"<<std::endl;
    auto fan_id = m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");
    auto dn2tmn = m_degenerated_mesh->getVariable<int,GMDS_NODE>("degMeshNode2TMeshNode");
    auto fan_angle = m_degenerated_mesh->getVariable<double,GMDS_NODE>("fan_angle");
    auto tj = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
    auto corner1 = m_degenerated_mesh->getVariable<int,GMDS_NODE>("corner");
    auto corner2 = m_t_mesh->getVariable<int,GMDS_NODE>("corner");
    auto constraints_nodes1 = m_degenerated_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto constraints_nodes2 = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto n2g = m_degenerated_mesh->getVariable<int,GMDS_NODE>("degenerated_node_to_group");
    auto iaut = m_t_mesh->getVariable<int,GMDS_FACE>("is_an_unwanted_triangle");
    auto iautv = m_t_mesh->getVariable<int,GMDS_NODE>("is_an_unwanted_triangle_vertex");
    auto of2tmf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("oldFace2tmeshFace");
    // Find the number of fans
    int nb_fans = -1;
    for (auto n_id:m_degenerated_mesh->nodes())
    {
        if (fan_id->value(n_id) > nb_fans)
            nb_fans = fan_id->value(n_id);
    }
    nb_fans += 1;
    // list of fans
    std::vector<std::vector<Face>> fans;
    for (int i = 0; i < nb_fans; i++)
    {
        std::vector<Face> tri_fan = fan(i);
        fans.push_back(tri_fan);
    }
    // Type of the fan (1 if there is 1 quad to put, 2 if 2 quads)
    std::vector<int> fan_type(fans.size());
    for (auto n_id:m_degenerated_mesh->nodes())
    {
        if (fan_id->value(n_id) >= 0)
        {
            if (fabs(fan_angle->value(n_id)-3.*M_PI/2.) < 0.1)
                fan_type[fan_id->value(n_id)] = 1;
            else
                fan_type[fan_id->value(n_id)] = 2;
        }
    }
    // Mark the faces belonging to a fan
    auto btaf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("belongs_to_a_fan");
    for (auto tri_fan:fans)
    {
        for (auto f:tri_fan)
        {
            btaf->set(f.id(),1);
        }
    }
    // Mark the fan vertices already added in the T-mesh
    std::vector<int> fan_vertex(nb_fans);
    for (int i = 0; i < nb_fans; i++)
        fan_vertex[i] = -1;

    // Build the nodes of the T-mesh
    for (auto n_id:m_degenerated_mesh->nodes())
    {
        int fanId = fan_id->value(n_id);
        if (fanId >= 0)
        {
            if (fan_vertex[fanId] == -1)
            {
                Node new_node = m_t_mesh->newNode(m_degenerated_mesh->get<Node>(n_id).point());
                fan_vertex[fanId] = new_node.id();
                dn2tmn->set(n_id,new_node.id());
                tj->set(new_node.id(),-1);
                corner2->set(new_node.id(),1);
                bool is_on_constraint = false;
                for (auto node:m_degenerated_nodes_groups[n2g->value(n_id)])
                {
                    if (constraints_nodes1->value(node.id()) == 1)
                    {
                        is_on_constraint = true;
                        break;
                    }
                }
                if (is_on_constraint)
                    constraints_nodes2->set(new_node.id(),1);
            }
            else
            {
                dn2tmn->set(n_id,fan_vertex[fanId]);
            }
        }
        else
        {
            Node n = m_degenerated_mesh->get<Node>(n_id);
            bool isInAFan = true;
            for (auto f:n.get<Face>())
            {
                if (btaf->value(f.id()) == 0)
                {
                    isInAFan = false;
                    break;
                }
            }
            if (!isInAFan)
            {
                Node new_node = m_t_mesh->newNode(m_degenerated_mesh->get<Node>(n_id).point());
                dn2tmn->set(n_id,new_node.id());
                tj->set(new_node.id(),-1);
                corner2->set(new_node.id(),corner1->value(n_id));
                constraints_nodes2->set(new_node.id(),constraints_nodes1->value(n_id));
            }
            else
            {
                dn2tmn->set(n_id,-1);
            }
        }
    }

    // Building the faces of the T-mesh
    for (auto f_id:m_degenerated_mesh->faces())
    {
        if (btaf->value(f_id) == 0)
        {
            std::vector<TCellID> nodes_id;
            Face f = m_degenerated_mesh->get<Face>(f_id);
            TCellID tri_ver;
            for (auto n:f.get<Node>())
            {
                bool already_seen = false;
                for (auto m_id:nodes_id)
                {
                    if (m_id == dn2tmn->value(n.id()))
                    {
                        already_seen = true;
                        tri_ver = m_id;
                        break;
                    }
                }
                if (!already_seen)
                    nodes_id.push_back(dn2tmn->value(n.id()));
            }
            Face new_face = m_t_mesh->newFace(nodes_id);
            if (nodes_id.size() == 3)
            {
                iaut->set(new_face.id(),1);
                iautv->set(tri_ver,1);
            }
            of2tmf->set(f_id,new_face.id());
        }
        else
        {
            of2tmf->set(f_id,-1);
        }
    }
    for (int i = 0; i < fans.size(); i++)
    {
        std::vector<Face> tri_fan = fans[i];
        std::vector<Node> outline_nodes = outline(tri_fan);
        if (fan_type[i] == 1)
        {
            std::vector<std::vector<Node>> sub_outlines = subOutlines1(outline_nodes,tri_fan,i);

            if (sub_outlines[1].size() < 2 || sub_outlines[2].size() < 2)
            {
                // It means that the fan is only one edge wide
                // Get the edge
                Edge fan_width;
                for (auto e1:sub_outlines[0][sub_outlines[0].size()-1].get<Edge>())
                {
                    for (auto e2:sub_outlines[3][0].get<Edge>())
                    {
                        if (e1.id() == e2.id())
                        {
                            fan_width = e1;
                            break;
                        }
                    }
                }
                // Get the face behind fan_width
                Face face_behind;
                for (auto f:fan_width.get<Face>())
                {
                    if (btaf->value(f.id()) == 0)
                    {
                        face_behind = f;
                        break;
                    }
                }
                // Create the new node
                Node middle_node = m_t_mesh->newNode((1./2.)*(fan_width.get<Node>()[0].point()+fan_width.get<Node>()[1].point()));
                // Create the face meshing the fan
                std::vector<TCellID> nodes_id;
                for (auto n:sub_outlines[0])
                {
                    nodes_id.push_back(dn2tmn->value(n.id()));
                }
                nodes_id.push_back(middle_node.id());
                for (int i = 0; i < sub_outlines[3].size()-1; i++)
                    nodes_id.push_back(dn2tmn->value(sub_outlines[3][i].id()));
                Face new_face = m_t_mesh->newFace(nodes_id);

                // Update the T-junctions
                for (int k = 1; k < sub_outlines[0].size()-1; k++)
                {
                    tj->set(dn2tmn->value(sub_outlines[0][k].id()),new_face.id());
                }
                for (int k = 1; k < sub_outlines[3].size()-1; k++)
                {
                    tj->set(dn2tmn->value(sub_outlines[3][k].id()),new_face.id());
                }

                // Update the face behind
                if (of2tmf->value(face_behind.id()) == -1)
                    throw GMDSException("BuildTMesh() : configuration not yet implemented");
                Face fb = m_t_mesh->get<Face>(of2tmf->value(face_behind.id()));
                nodes_id.clear();
                bool added = false;
                std::vector<Node> adj_nodes = fb.get<Node>();
                for (int i = 0; i < adj_nodes.size(); i++)
                {
                    if (added)
                        nodes_id.push_back(adj_nodes[i].id());
                    else
                    {
                        if (adj_nodes[i].id() == dn2tmn->value(fan_width.get<Node>()[0].id()) && adj_nodes[(i+1)%(adj_nodes.size())].id() == dn2tmn->value(fan_width.get<Node>()[1].id()))
                        {
                            nodes_id.push_back(adj_nodes[i].id());
                            nodes_id.push_back(middle_node.id());
                            added = true;
                        }
                        else if (adj_nodes[i].id() == dn2tmn->value(fan_width.get<Node>()[1].id()) && adj_nodes[(i+1)%(adj_nodes.size())].id() == dn2tmn->value(fan_width.get<Node>()[0].id()))
                        {
                            nodes_id.push_back(adj_nodes[i].id());
                            nodes_id.push_back(middle_node.id());
                            added = true;
                        }
                        else
                            nodes_id.push_back(adj_nodes[i].id());
                    }
                    
                }
                new_face = m_t_mesh->newFace(nodes_id);
                m_t_mesh->deleteFace(fb);
                tj->set(middle_node.id(),new_face.id());
            }

            else
            {
                std::vector<TCellID> nodes_id;
                for (auto n:outline_nodes)
                {
                    nodes_id.push_back(dn2tmn->value(n.id()));
                }
                Face new_face = m_t_mesh->newFace(nodes_id);

                // Update the T-junctions
                for (auto suboutline:sub_outlines)
                {
                    for (int k = 1; k < suboutline.size()-1; k++)
                    {
                        tj->set(dn2tmn->value(suboutline[k].id()),new_face.id());
                    }
                }
            }

        }
        if (fan_type[i] == 2)
        {
            std::vector<std::vector<Node>> sub_outlines = subOutlines2(outline_nodes,tri_fan,i);
            
            // for (auto sub:sub_outlines)
            // {
            //     std::cout<<"TEST "<<std::endl;
            //     for (auto n:sub)
            //         std::cout<<"test "<<n.id()<<std::endl;
            // }

            
            std::vector<TCellID> nodes_id;
            for (int j = 0; j < sub_outlines[0].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[0][j].id()));
            }
            for (int j = 1; j < sub_outlines[1].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[1][j].id()));
            }
            for (int j = 1; j < sub_outlines[2].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[2][j].id()));
            }
            Face new_face = m_t_mesh->newFace(nodes_id);

            // Update the T-junctions
            for (int j = 0; j <= 2; j++)
            {
                for (int k = 1; k < sub_outlines[j].size()-1; k++)
                {
                    tj->set(dn2tmn->value(sub_outlines[j][k].id()),new_face.id());
                }
            }


            nodes_id.clear();
            for (int j = 0; j < sub_outlines[3].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[3][j].id()));
            }
            for (int j = 1; j < sub_outlines[4].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[4][j].id()));
            }
            for (int j = 1; j < sub_outlines[5].size(); j++)
            {
                nodes_id.push_back(dn2tmn->value(sub_outlines[5][j].id()));
            }
            new_face = m_t_mesh->newFace(nodes_id);

            // Update the T-junctions
            for (int j = 3; j <= 5; j++)
            {
                for (int k = 1; k < sub_outlines[j].size()-1; k++)
                {
                    tj->set(dn2tmn->value(sub_outlines[j][k].id()),new_face.id());
                }
            }
        }
    }
    m_degenerated_mesh->deleteVariable(GMDS_FACE,btaf);
    m_degenerated_mesh->deleteVariable(GMDS_FACE,of2tmf);
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::buildDegeneratedNodesGroups()
{
    std::cout<<"> Building groups of degenerated nodes"<<std::endl;
    auto fan_id = m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");
    auto n2g = m_degenerated_mesh->getVariable<int,GMDS_NODE>("degenerated_node_to_group");
    for (auto n_id:m_degenerated_mesh->nodes())
        n2g->set(n_id,-1);
    std::vector<std::vector<Node>> groups;
    int group_id = 0;
    for (auto n_id:m_degenerated_mesh->nodes())
    {
        std::vector<Node> group;
        if (fan_id->value(n_id) == -1)
        {
            Node n = m_degenerated_mesh->get<Node>(n_id);
            group.push_back(n);
            groups.push_back(group);
            n2g->set(n_id,group_id);
            group_id += 1;
        }
        else
        {
            if (n2g->value(n_id) < 0)
            {
                int fanId = fan_id->value(n_id);
                for (auto m_id:m_degenerated_mesh->nodes())
                {
                    if (fan_id->value(m_id) == fanId)
                    {
                        Node n = m_degenerated_mesh->get<Node>(m_id);
                        group.push_back(n);
                        n2g->set(m_id,group_id);
                    }
                }
                groups.push_back(group);
                group_id += 1;
            }
        }
    }

    m_degenerated_nodes_groups = groups;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> TrianglesRemover::subOutlines1(std::vector<Node> AOutline, std::vector<Face> AFan, int AFanId)
{
    auto fan_id = m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");
    auto bttf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("belongs_to_the_fan");
    auto separatrix = m_degenerated_mesh->getVariable<int,GMDS_EDGE>("belongs_to_a_separatrix");
    for (auto f:AFan)
    {
        bttf->set(f.id(),1);
    }
    // Find the triangle vertex
    int pos;
    for (int i = 0; i < AOutline.size(); i++)
    {
        if (fan_id->value(AOutline[i].id()) == AFanId)
        {
            pos = i;
            break;
        }
    }
    std::vector<std::vector<Node>> sub_outlines;
    std::vector<Node> suboutline,suboutline2;
    int i = pos;
    bool Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        int NbFanFaces = 0;
        for (auto f:AOutline[i].get<Face>())
        {
            if (bttf->value(f.id()) == 1)
            {
                NbFanFaces += 1;
            }
        }
        if (NbFanFaces == 1)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    sub_outlines.push_back(suboutline);

    suboutline.clear();
    Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        int NbFanFaces = 0;
        for (auto f:AOutline[i].get<Face>())
        {
            if (bttf->value(f.id()) == 1)
            {
                NbFanFaces += 1;
            }
        }
        if (NbFanFaces == 1)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    int half = suboutline.size()/2;

    // Look for singularities
    bool singu_found = false;
    int singu_pos;
    for (int j = 1; j < suboutline.size()-1; j++)
    {
        int NbNonSingu = 0;
        for (auto e:suboutline[j].get<Edge>())
        {
            if (separatrix->value(e.id()) == 0)
            {
                NbNonSingu += 1;
            }
        }
        if (NbNonSingu == 0)
        {
            singu_found = true;
            singu_pos = j;
            break;
        }
    }


    if (singu_found)
    {
        for (int j = 0; j <= singu_pos; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        suboutline2.clear();
        for (int j = singu_pos; j < suboutline.size(); j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }
    else
    {
        for (int j = 0; j <= half; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        suboutline2.clear();
        for (int j = half; j < suboutline.size(); j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }


    suboutline.clear();
    Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        if (i == pos)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    sub_outlines.push_back(suboutline);

    m_degenerated_mesh->deleteVariable(GMDS_FACE,bttf);
    return sub_outlines;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> TrianglesRemover::subOutlines2(std::vector<Node> AOutline, std::vector<Face> AFan, int AFanId)
{
    auto fan_id = m_degenerated_mesh->getVariable<int,GMDS_NODE>("fan_id");
    auto bttf = m_degenerated_mesh->newVariable<int,GMDS_FACE>("belongs_to_the_fan");
    auto separatrix = m_degenerated_mesh->getVariable<int,GMDS_EDGE>("belongs_to_a_separatrix");
    for (auto f:AFan)
    {
        bttf->set(f.id(),1);
    }
    // Find the triangle vertex
    int pos;
    for (int i = 0; i < AOutline.size(); i++)
    {
        if (fan_id->value(AOutline[i].id()) == AFanId)
        {
            pos = i;
            break;
        }
    }
    std::vector<std::vector<Node>> sub_outlines;
    std::vector<Node> suboutline,suboutline2;
    int i = pos;
    bool Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        int NbFanFaces = 0;
        for (auto f:AOutline[i].get<Face>())
        {
            if (bttf->value(f.id()) == 1)
            {
                NbFanFaces += 1;
            }
        }
        if (NbFanFaces == 1)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    sub_outlines.push_back(suboutline);

    suboutline.clear();
    Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        int NbFanFaces = 0;
        for (auto f:AOutline[i].get<Face>())
        {
            if (bttf->value(f.id()) == 1)
            {
                NbFanFaces += 1;
            }
        }
        if (NbFanFaces == 1)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    int quarter = suboutline.size()/4;
    int half = suboutline.size()/2;

    // Look for singularities
    bool singu_found = false;
    int singu_pos;
    for (int j = 1; j < half; j++)
    {
        int NbNonSingu = 0;
        for (auto e:suboutline[j].get<Edge>())
        {
            if (separatrix->value(e.id()) == 0)
            {
                NbNonSingu += 1;
            }
        }
        if (NbNonSingu == 0)
        {
            singu_found = true;
            singu_pos = j;
            break;
        }
    }

    if (singu_found)
    {
        for (int j = 0; j <= singu_pos; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        
        suboutline2.clear();
        for (int j = singu_pos; j <= half; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }
    else
    {
        for (int j = 0; j <= quarter; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        
        suboutline2.clear();
        for (int j = quarter; j <= half; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }


    // Look for singularities
    singu_found = false;
    for (int j = half+1; j < suboutline.size()-1; j++)
    {
        int NbNonSingu = 0;
        for (auto e:suboutline[j].get<Edge>())
        {
            if (separatrix->value(e.id()) == 0)
            {
                NbNonSingu += 1;
            }
        }
        if (NbNonSingu == 0)
        {
            singu_found = true;
            singu_pos = j;
            break;
        }
    }

    if (singu_found)
    {
        suboutline2.clear();
        for (int j = half; j <= singu_pos; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        
        suboutline2.clear();
        for (int j = singu_pos; j < suboutline.size(); j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }
    else
    {
        suboutline2.clear();
        for (int j = half; j <= half+quarter; j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
        
        suboutline2.clear();
        for (int j = half+quarter; j < suboutline.size(); j++)
        {
            suboutline2.push_back(suboutline[j]);
        }
        sub_outlines.push_back(suboutline2);
    }


    suboutline.clear();
    Continue = true;
    while (Continue)
    {
        suboutline.push_back(AOutline[i]);
        i = i+1;
        if (i == AOutline.size())
            i = 0;
        if (i == pos)
        {
            Continue = false;
        }
    }
    suboutline.push_back(AOutline[i]);
    sub_outlines.push_back(suboutline);

    m_degenerated_mesh->deleteVariable(GMDS_FACE,bttf);
    return sub_outlines;
}

/*----------------------------------------------------------------------------*/
Mesh TrianglesRemover::getTMesh()
{
	return *m_t_mesh;
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::setTMeshConnectivity()
{
	std::cout<<"> Setting T-mesh connectivity"<<std::endl;
	MeshDoctor doc(m_t_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::markInternalConstraintsOnEdges()
{
	std::cout<<"> Marking internal constraints on edges of the T-mesh"<<std::endl;
    auto constraints_nodes1 = m_degenerated_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto constraints_nodes2 = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto constraints_edges1 = m_degenerated_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto constraints_edges2 = m_t_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto dmn2tmn = m_degenerated_mesh->getVariable<int,GMDS_NODE>("degMeshNode2TMeshNode");
    for (auto e_id:m_degenerated_mesh->edges())
    {
        if (constraints_edges1->value(e_id) == 1)
        {
            Edge e = m_degenerated_mesh->get<Edge>(e_id);
            Node n1 = e.get<Node>()[0];
            Node n2 = e.get<Node>()[1];
            if (dmn2tmn->value(n1.id()) >= 0 && dmn2tmn->value(n2.id()) >= 0 && dmn2tmn->value(n1.id()) != dmn2tmn->value(n2.id()))
            {
                Node m1 = m_t_mesh->get<Node>(dmn2tmn->value(n1.id()));
                Node m2 = m_t_mesh->get<Node>(dmn2tmn->value(n2.id()));
                try
                {
                    Edge e2 = getEdge(m1,m2);
                    constraints_edges2->set(e2.id(),1);
                }
                catch (GMDSException& e){}
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::transformUnwantedTrianglesIntoQuads()
{
    std::cout<<"> Transforming the unwanted triangles into quads"<<std::endl;
    auto iaut = m_t_mesh->getVariable<int,GMDS_FACE>("is_an_unwanted_triangle");
    auto iautv = m_t_mesh->getVariable<int,GMDS_NODE>("is_an_unwanted_triangle_vertex");
    auto constraints_nodes = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto constraints_edges = m_t_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto tj = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
    auto iann = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_new_node");
    // Mark edges belonging to triangles and triangles vertices
    auto tri_edge = m_t_mesh->newVariable<int,GMDS_EDGE>("is_a_tri_edge");
    for (auto n_id:m_t_mesh->nodes())
    {
        if (iautv->value(n_id) == 1)
        {
            Node vertex = m_t_mesh->get<Node>(n_id);
            std::vector<Edge> adj_edges = vertex.get<Edge>();
            adj_edges = sortEdges(vertex,adj_edges);
            for (auto e:adj_edges)
            {
                bool is_a_tri_edge = false;
                for (auto f:e.get<Face>())
                {
                    if (iaut->value(f.id()) == 1)
                    {
                        is_a_tri_edge = true; 
                        break;
                    }
                }
                if (is_a_tri_edge)
                    tri_edge->set(e.id(),1);
            }
            // Identify the fan of tri edges
            std::vector<Edge> edges_fan;
            int start_pos;
            for (int i = 0; i < adj_edges.size(); i++)
            {
                if (tri_edge->value(adj_edges[i].id()) == 0 && tri_edge->value(adj_edges[(i+1)%(adj_edges.size())].id()) == 1 && tri_edge->value(adj_edges[(i+2)%(adj_edges.size())].id()) == 1)
                {
                    start_pos = i;
                    break;
                }
            }
            edges_fan.push_back(adj_edges[start_pos]);
            int incr = (start_pos+1)%(adj_edges.size());
            while (tri_edge->value(adj_edges[incr].id()) == 1)
            {
                edges_fan.push_back(adj_edges[incr]);
                incr = (incr + 1)%(adj_edges.size());
            }
            edges_fan.push_back(adj_edges[incr]);
            // Check that the fan is well oriented
            if (constraints_edges->value(edges_fan[1].id()) == 1)
                std::reverse(edges_fan.begin(),edges_fan.end());

            // Create the new points and quads
            math::Vector dir = edge2vec(edges_fan[0],vertex);
            std::vector<TCellID> new_nodes;
            Node new_node;
            math::Point P;
            for (int i = 1; i < edges_fan.size()-2; i++)
            {
                P = vertex.point() + (double(i)/double(edges_fan.size()-2))*dir;
                new_node = m_t_mesh->newNode(P);
                new_nodes.push_back(new_node.id());
            }

            std::vector<TCellID> above_nodes;
            Face neighbouring_face = getCommonFace(edges_fan[0],edges_fan[1]);
            TCellID right_node,left_node,last_node_added;
            if (edges_fan[0].get<Node>()[0].id() == vertex.id())
                right_node = edges_fan[0].get<Node>()[1].id();
            else
                right_node = edges_fan[0].get<Node>()[0].id();
            if (edges_fan[1].get<Node>()[0].id() == vertex.id())
                left_node = edges_fan[1].get<Node>()[1].id();
            else
                left_node = edges_fan[1].get<Node>()[0].id();
            for (auto n:neighbouring_face.get<Node>())
            {
                if (n.id() != vertex.id() && n.id() != right_node && n.id() != left_node)
                {
                    above_nodes.push_back(n.id());
                    break;
                }
            }
            above_nodes.push_back(left_node);
            last_node_added = left_node;
            for (int i = 1; i < edges_fan.size()-2; i++)
            {
                Face f = getCommonFace(edges_fan[i],edges_fan[i+1]);
                for (auto n:f.get<Node>())
                {
                    if (n.id() != vertex.id() && n.id() != last_node_added)
                    {
                        above_nodes.push_back(n.id());
                        last_node_added = n.id();
                        break;
                    }
                }
                m_t_mesh->deleteFace(f);
            }
            m_t_mesh->deleteFace(neighbouring_face);
            std::vector<TCellID> below_nodes;
            below_nodes.push_back(vertex.id());
            for (auto m_id:new_nodes)
                below_nodes.push_back(m_id);
            below_nodes.push_back(right_node);
            // Create the new quads
            Face new_face;
            std::vector<TCellID> new_face_ids;
            for (int i = 0; i < below_nodes.size()-1; i++)
            {
                new_face_ids = {below_nodes[i],below_nodes[i+1],above_nodes[above_nodes.size()-2-i],above_nodes[above_nodes.size()-1-i]};
                new_face = m_t_mesh->newFace(new_face_ids);
            }
            // Change the behind face
            Face behind_face;
            for (auto f:edges_fan[0].get<Face>())
            {
                if (f.id() != neighbouring_face.id())
                {
                    behind_face = f;
                    break;
                }
            }
            new_face_ids.clear();
            std::vector<Node> adj_nodes = behind_face.get<Node>();
            bool already_added = false;
            for (int i = 0; i < adj_nodes.size(); i++)
            {
                if (already_added)
                    new_face_ids.push_back(adj_nodes[i].id());
                else
                {
                    if (adj_nodes[i].id() == vertex.id() && adj_nodes[(i+1)%adj_nodes.size()].id() == right_node)
                    {
                        new_face_ids.push_back(adj_nodes[i].id());
                        for (auto id:new_nodes)
                            new_face_ids.push_back(id);
                        already_added = true;
                    }
                    else if (adj_nodes[(i+1)%adj_nodes.size()].id() == vertex.id() && adj_nodes[i].id() == right_node)
                    {
                        new_face_ids.push_back(adj_nodes[i].id());
                        for (int j = 0; j < new_nodes.size(); j++)
                            new_face_ids.push_back(new_nodes[new_nodes.size()-1-j]);
                        already_added = true;
                    }
                    else
                        new_face_ids.push_back(adj_nodes[i].id());
                }
            }
            new_face = m_t_mesh->newFace(new_face_ids);
            for (auto id:new_nodes)
                tj->set(id,new_face.id());
            m_t_mesh->deleteFace(behind_face);
            if (constraints_edges->value(edges_fan[0].id()) == 1)
            {
                for (auto id:new_nodes)
                {
                    constraints_nodes->set(id,1);
                    iann->set(id,1);
                }                
            }
        }
    }


    m_t_mesh->deleteVariable(GMDS_EDGE,tri_edge);
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::buildFinalTMesh()
{
	std::cout<<"> Building the final T-mesh"<<std::endl;
    auto tj1 = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
    auto c1 = m_t_mesh->getVariable<int,GMDS_NODE>("corner");
    auto btic1 = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
	auto tj2 = m_final_t_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
    auto c2 = m_final_t_mesh->getVariable<int,GMDS_NODE>("corner");
    auto btic2 = m_final_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto of2nf = m_t_mesh->newVariable<int,GMDS_FACE>("oldFace2newFace");
    auto iann1 = m_t_mesh->getVariable<int,GMDS_NODE>("is_a_new_node");
    auto iann2 = m_final_t_mesh->getVariable<int,GMDS_NODE>("is_a_new_node");
    // Build nodes
    Node n,new_node;
    for (auto n_id:m_t_mesh->nodes())
    {
        n = m_t_mesh->get<Node>(n_id);
        new_node = m_final_t_mesh->newNode(n.point());
        c2->set(new_node.id(),c1->value(n_id));
        btic2->set(new_node.id(),btic1->value(n_id));
        iann2->set(new_node.id(),iann1->value(n_id));
    }
    // Build faces
    Face f,new_face;
    for (auto f_id:m_t_mesh->faces())
    {
        f = m_t_mesh->get<Face>(f_id);
        std::vector<TCellID> ids;
        for (auto n:f.get<Node>())
            ids.push_back(n.id());
        new_face = m_final_t_mesh->newFace(ids);
        of2nf->set(f_id,new_face.id());
    }
    // Update the T-junctions
    for (auto n_id:m_final_t_mesh->nodes())
    {
        if (tj1->value(n_id) == -1)
            tj2->set(n_id,-1);
        else
            tj2->set(n_id,of2nf->value(tj1->value(n_id)));
    }

    m_t_mesh->deleteVariable(GMDS_FACE,of2nf);
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::writeFinalTMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the final T-mesh"<<std::endl;
	IGMeshIOService ioService(m_final_t_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
Mesh TrianglesRemover::getFinalTMesh()
{
	return *m_final_t_mesh;
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::setFinalTMeshConnectivity()
{
	std::cout<<"> Setting the final T-mesh connectivity"<<std::endl;
	MeshDoctor doc(m_final_t_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void TrianglesRemover::markInternalConstraintsOnFinalEdges()
{
	std::cout<<"> Marking internal constraints on the edges of the final T-mesh"<<std::endl;
	auto ic1 = m_t_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto btic1 = m_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto ic2 = m_final_t_mesh->getVariable<int,GMDS_EDGE>("internal_constraint");
    auto btic2 = m_final_t_mesh->getVariable<int,GMDS_NODE>("belong_to_an_internal_constraint");
    auto iann = m_final_t_mesh->getVariable<int,GMDS_NODE>("is_a_new_node");
    for (auto e_id:m_final_t_mesh->edges())
    {
        Edge e1 = m_final_t_mesh->get<Edge>(e_id);
        Node n1 = m_t_mesh->get<Node>(e1.get<Node>()[0].id());
        Node n2 = m_t_mesh->get<Node>(e1.get<Node>()[1].id());
        try
        {
            Edge e2 = getEdge(n1,n2);
            ic2->set(e1.id(),ic1->value(e2.id()));
        }
        catch (GMDSException& e){}
        if (iann->value(e1.get<Node>()[0].id()) == 1 || iann->value(e1.get<Node>()[1].id()) == 1)
        {
            if (btic2->value(e1.get<Node>()[0].id()) == 1 && btic2->value(e1.get<Node>()[1].id()) == 1)
                ic2->set(e1.id(),1);
        }

    }
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/