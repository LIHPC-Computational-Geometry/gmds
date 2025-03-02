#include "gmds/medialaxis/MinDelaunayCleaner.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
MinDelaunayCleaner::MinDelaunayCleaner(gmds::Mesh &AMesh)
{
	// Input minimal Delaunay trianglulation
	m_mesh = &AMesh;
    // Face type
    m_mesh->newVariable<int,GMDS_FACE>("face_type");
    // Medial points
    m_mesh->newVariable<math::Point,GMDS_FACE>("medial_point");
    // Medial radii
    m_mesh->newVariable<double,GMDS_FACE>("medial_radius");
    // Medial radii
    m_mesh->newVariable<int,GMDS_FACE>("medial_branch_id");
    // Faces to delete are marked with 1
    m_mesh->newVariable<int,GMDS_FACE>("is_to_keep");
    // Small edges are marked with 1
    m_mesh->newVariable<int,GMDS_EDGE>("is_small");

	// Cleaned mesh
	m_cleaned_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));

}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::setFacesTypes()
{
    auto face_type = m_mesh->getVariable<int,GMDS_FACE>("face_type");
    for (auto f_id:m_mesh->faces())
    {
        Face f = m_mesh->get<Face>(f_id);
        int type = 0;
        for (auto e:f.get<Edge>())
            type += e.get<Face>().size()-1;
        face_type->set(f_id,type);
    }
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::computeMedialPointsAndRadii()
{
    auto mp = m_mesh->getVariable<math::Point,GMDS_FACE>("medial_point");
    auto mr = m_mesh->getVariable<double,GMDS_FACE>("medial_radius");
    for (auto f_id:m_mesh->faces())
    {
        Face f = m_mesh->get<Face>(f_id);
        math::Triangle T(f.get<Node>()[0].point(),f.get<Node>()[1].point(),f.get<Node>()[2].point());
        math::Point MP = T.getCircumcenter();
        double r = MP.distance(f.get<Node>()[0].point());
        mp->set(f_id,MP);
        mr->set(f_id,r);
    }
}

/*----------------------------------------------------------------------------*/
std::vector<Face> MinDelaunayCleaner::medialBranch(Face &AF)
{
    auto face_type = m_mesh->getVariable<int,GMDS_FACE>("face_type");
    std::vector<Face> medial_branch;
    if (face_type->value(AF.id()) == 3)
        return medial_branch;
    std::queue<Face> toAdd;
    Face current = AF;
    toAdd.push(current);
    auto visited = m_mesh->newVariable<int,GMDS_FACE>("visited");
    visited->set(current.id(),1);
    while(!toAdd.empty())
    {
        current = toAdd.front();
        toAdd.pop();
        medial_branch.push_back(current);
        for (auto e:current.get<Edge>())
        {
            for (auto f:e.get<Face>())
            {
                if (f.id() != current.id() && face_type->value(f.id()) != 3 && visited->value(f.id()) == 0)
                {
                    toAdd.push(f);
                    visited->set(f.id(),1);
                }
            }
        }
    }
    m_mesh->deleteVariable(GMDS_FACE,visited);
    return medial_branch;
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::setBranchID()
{
    auto id = m_mesh->getVariable<int,GMDS_FACE>("medial_branch_id");
    for (auto f_id:m_mesh->faces())
        id->set(f_id,-1);
    int ID = 0;
    for (auto f_id:m_mesh->faces())
    {
        if (id->value(f_id) == -1)
        {
            Face f = m_mesh->get<Face>(f_id);
            std::vector<Face> branch = medialBranch(f);
            for (auto f1:branch)
                id->set(f1.id(),ID);
            ID += 1;
        }
    }
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::markSmallEdges()
{
    auto is_small = m_mesh->getVariable<int,GMDS_EDGE>("is_small");
    double mean_len = 0.;
    for (auto e_id:m_mesh->edges())
    {
        Edge e = m_mesh->get<Edge>(e_id);
        mean_len += e.length();
    }
    mean_len = mean_len/double(m_mesh->getNbEdges());
    for (auto e_id:m_mesh->edges())
    {
        Edge e = m_mesh->get<Edge>(e_id);
        // if (e.length() < mean_len/2.)
        if (e.length() < mean_len/1.25)
            is_small->set(e_id,1);
    }
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::markFacesToDelete()
{
    auto is_to_keep = m_mesh->getVariable<int,GMDS_FACE>("is_to_keep");
    auto face_type = m_mesh->getVariable<int,GMDS_FACE>("face_type");
    auto is_small = m_mesh->getVariable<int,GMDS_EDGE>("is_small");
    // Mark big intersection triangles
    auto big_itri = m_mesh->newVariable<int,GMDS_FACE>("is_a_big_itri");
    for (auto f_id:m_mesh->faces())
    {
        if (face_type->value(f_id) == 3)
        {
            Face f = m_mesh->get<Face>(f_id);
            bool isBig = true;
            for (auto e:f.get<Edge>())
            {
                if (is_small->value(e.id()) == 1)
                {
                    isBig = false;
                    break;
                }
            }
            if (isBig)
                big_itri->set(f_id,1);
        }
    }
    // Check if there is at least one big intersection triangle. If not, artificially add one
    bool valid = false;
    for (auto f_id:m_mesh->faces())
    {
        if (big_itri->value(f_id) == 1)
        {
            valid = true;
            break;
        }
    }
    if (!valid)
    {
        for (auto f_id:m_mesh->faces())
        {
            big_itri->set(f_id,1);
            break;
        }
    }
    // Find the triangles to keep from the big intersection triangles
    for (auto f_id:m_mesh->faces())
    {
        if (big_itri->value(f_id) == 1)
        {
            Face f = m_mesh->get<Face>(f_id);
            is_to_keep->set(f_id,1);
            for (auto e:f.get<Edge>())
            {
                for (auto f1:e.get<Face>())
                {
                    if (f1.id() != f.id() && is_to_keep->value(f1.id()) == 0)
                    {
                        Face current_face = f1;
                        Edge current_edge = e;
                        is_to_keep->set(f1.id(),1);
                        while (big_itri->value(current_face.id()) == 0 && face_type->value(current_face.id()) != 1)
                        {
                            // We go through the biggest edge of current_face
                            Edge biggest_edge;
                            double l = 0;
                            for (auto e1:current_face.get<Edge>())
                            {
                                if (e1.id() != current_edge.id() && e1.length() > l && e1.get<Face>().size() == 2)
                                {
                                    biggest_edge = e1;
                                    l = e1.length();
                                }
                            }
                            current_edge = biggest_edge;
                            for (auto f2:current_edge.get<Face>())
                            {
                                if (f2.id() != current_face.id())
                                {
                                    current_face = f2;
                                    break;
                                }
                            }
                            is_to_keep->set(current_face.id(),1);
                        }
                    }
                }
            }
        }
    }


    m_mesh->deleteVariable(GMDS_FACE,big_itri);
    
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::buildCleanedMesh()
{
    auto is_to_keep = m_mesh->getVariable<int,GMDS_FACE>("is_to_keep");
    // Mark the nodes to keep
    auto node_to_keep = m_mesh->newVariable<int,GMDS_NODE>("node_to_keep");
    for (auto f_id:m_mesh->faces())
    {
        if (is_to_keep->value(f_id) == 1)
        {
            Face f = m_mesh->get<Face>(f_id);
            for (auto n:f.get<Node>())
                node_to_keep->set(n.id(),1);
        }
    }

    // Build the new nodes
    auto old2new = m_mesh->newVariable<TCellID,GMDS_NODE>("old2new");
    for (auto n_id:m_mesh->nodes())
    {
        if (node_to_keep->value(n_id) == 1)
        {
            Node n = m_mesh->get<Node>(n_id);
            Node newNode = m_cleaned_mesh->newNode(n.point());
            old2new->set(n_id,newNode.id());
        }
    }

    // Build the new faces
    for (auto f_id:m_mesh->faces())
    {
        if (is_to_keep->value(f_id) == 1)
        {
            Face f = m_mesh->get<Face>(f_id);
            std::vector<TCellID> nodes;
            for (auto n:f.get<Node>())
                nodes.push_back(old2new->value(n.id()));
            m_cleaned_mesh->newFace(nodes);
        }
    }
    m_mesh->deleteVariable(GMDS_NODE,node_to_keep);
    m_mesh->deleteVariable(GMDS_NODE,old2new);
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::setCleanedMeshConnectivity()
{
    MeshDoctor doc(m_cleaned_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
Mesh MinDelaunayCleaner::getCleanedMesh()
{
    return *m_cleaned_mesh;
}

/*----------------------------------------------------------------------------*/
void MinDelaunayCleaner::writeCleanedMesh(std::string AFileName)
{
    IGMeshIOService ioService(m_cleaned_mesh);
    VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/