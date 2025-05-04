/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Face.h>
#include <gmds/ig/FaceContainer.h>
#include <gmds/ig/Mesh.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Quadrilateral.h>
#include <gmds/math/Segment.h>
#include <gmds/utils/Exception.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
    Face::Face()
            : Cell(nullptr, GMDS_FACE, NullID)
            , m_faces_container(nullptr)
            , m_type_id(NullID)
    {
    }
/*----------------------------------------------------------------------------*/
    TInt
    Face::nbNodes() const
    {
            TInt nb = 0;
            m_faces_container->getNodesData(m_id, nb);
            return nb;
    }
/*----------------------------------------------------------------------------*/
    TInt
    Face::nbEdges() const
    {
            TInt nb = 0;
            m_faces_container->getEdgesData(m_id, nb);
            return nb;
    }
/*----------------------------------------------------------------------------*/
    TInt
    Face::nbFaces() const
    {
            TInt nb = 0;
            m_faces_container->getFacesData(m_id, nb);
            return nb;
    }
/*----------------------------------------------------------------------------*/
    TInt
    Face::nbRegions() const
    {
            TInt nb = 0;
            m_faces_container->getRegionsData(m_id, nb);
            return nb;
    }
/*----------------------------------------------------------------------------*/
    math::Vector3d
    Face::normal() const
    {
            math::Vector3d n;
            std::vector<Node> nodes = this->get<Node>();
            auto nb_nodes = nodes.size();

            if (nb_nodes == 3) {
                    math::Point p1 = nodes[0].point();
                    math::Point p2 = nodes[1].point();
                    math::Point p3 = nodes[2].point();

                    math::Vector3d v1=p2-p1;
                    math::Vector3d v3=p3-p1;

                    n = (v1.cross(v3));
                    n.normalize();
            } else if (nb_nodes > 3){
                    math::Vector3d n_temp({0,0,0});
		              int tmp_nb=0;
                    for (auto i = 0; i < nb_nodes; i++) {
			             		 math::Point p1 = nodes[(i==0)?nb_nodes-1:i - 1].point();
                            math::Point p2 = nodes[i].point();
                            math::Point p3 = nodes[(i + 1) % nb_nodes].point();
			                   math::Vector3d v1=p2-p1;
			                   math::Vector3d v3=p3-p1;

			                   math::Vector3d ni = (v1.cross(v3));
			                   n_temp+=ni;
			                   tmp_nb++;
                    }
                    n_temp /= tmp_nb;
                    n_temp.normalize();
                    n = n_temp;
            }
            return n;
    }
/*----------------------------------------------------------------------------*/
    math::Vector3d
    Face::normal(const Node& ANode) const
    {
            math::Vector3d n;
            std::vector<Node> nodes = this->get<Node>();

            // first find ANode among the faces nodes
            unsigned int nodeIndex = NullID;

            for (unsigned int iNode = 0; iNode < nodes.size(); iNode++) {
                    if (ANode == nodes[iNode]) {
                            nodeIndex = iNode;
                    }
            }

            if (nodeIndex == NullID) {
                    throw GMDSException("math::Vector Face::normal could not find ANode.");
            }

            Node n1 = nodes[nodeIndex];
            Node n2 = nodes[(nodeIndex + 1) % nodes.size()];
            Node n3 = nodes[(nodeIndex + 2) % nodes.size()];

            math::Vector3d v1({n2.X() - n1.X(), n2.Y() - n1.Y(), n2.Z() - n1.Z()});
            math::Vector3d v3({n3.X() - n1.X(), n3.Y() - n1.Y(), n3.Z() - n1.Z()});

            n = (v1.cross(v3));
            n.normalize();
            return n;
    }
/*----------------------------------------------------------------------------*/
    TCoord
    Face::area() const {
            math::Vector3d n;
            std::vector<Node> nodes = this->get<Node>();
            auto nb_nodes = nodes.size();

            if (nb_nodes == 3) {

                    Node n1 = nodes[0];
                    Node n2 = nodes[1];
                    Node n3 = nodes[2];

                    math::Vector3d v1({n2.X() - n1.X(), n2.Y() - n1.Y(), n2.Z() - n1.Z()});
                    math::Vector3d v3({n3.X() - n1.X(), n3.Y() - n1.Y(), n3.Z() - n1.Z()});

                    n = (v1.cross(v3));

                    return 0.5 * n.norm();
            } else if (nb_nodes == 4) {
                    math::Point p1 = nodes[0].point();
                    math::Point p2 = nodes[1].point();
                    math::Point p3 = nodes[2].point();
                    math::Point p4 = nodes[3].point();
                    math::Point c = center();
                    math::Triangle t1(p1,p2,c);
                    math::Triangle t2(p2,p3,c);
                    math::Triangle t3(p3,p4,c);
                    math::Triangle t4(p4,p1,c);

                    return t1.area()+t2.area()+t3.area()+t4.area();
            } else {
                    throw GMDSException("Not yet implemented!");

            }
    }
/*----------------------------------------------------------------------------*/
    math::Point
    Face::center() const
    {
            TCoord p_coords[3] = {0.0, 0.0, 0.0};

            std::vector<Node> nodes = this->get<Node>();
            auto nb_nodes = nodes.size();

            for (auto i = 0; i < nb_nodes; i++) {
                    Node n = nodes[i];
                    p_coords[0] += n.X();
                    p_coords[1] += n.Y();
                    p_coords[2] += n.Z();
            }

            p_coords[0] = p_coords[0] / nodes.size();
            p_coords[1] = p_coords[1] / nodes.size();
            p_coords[2] = p_coords[2] / nodes.size();

            math::Point p(p_coords[0], p_coords[1], p_coords[2]);

            return p;
    }
/*----------------------------------------------------------------------------*/
    double
    Face::computeScaledJacobian2D() const
    {
            std::vector<Node> nodes = this->get<Node>();

            switch (this->type()) {
                    case GMDS_TRIANGLE:
                            throw GMDSException("Region::computeScaledJacobian2D not implemented yet for triangles.");
                    break;
                    case GMDS_QUAD: {
                            gmds::math::Quadrilateral quad(nodes[0].point(), nodes[1].point(), nodes[2].point(),
                                                           nodes[3].point());
                            return quad.computeScaledJacobian2D();
                    } break;
                    default:
                            throw GMDSException("Region::computeScaledJacobian2D not implemented yet for this value type.");
                    break;
            }
        throw GMDSException("Region::computeScaledJacobian2D not implemented yet for this value type.");
    }
/*----------------------------------------------------------------------------*/
    math::Point
    Face::project(const math::Point& AP) const
    {
            if (m_type != GMDS_TRIANGLE && m_type != GMDS_QUAD)
                    throw GMDSException("Face::project only implemented for triangular and quadrangular faces");

            std::vector<math::Triangle> ts;
            if (m_type == GMDS_TRIANGLE){
                std::vector<Node> nodes = this->get<Node>();
                math::Point p0 = nodes[0].point();
                math::Point p1 = nodes[1].point();
                math::Point p2 = nodes[2].point();
                // we value the projected point
                math::Point X = math::Plane(p0, p1, p2).project(AP);

                TCoord x, y, z;
                math::Point::computeBarycentric(p0, p1, p2, X, x, y, z);

                if (x < 0.0) {
                    if (y < 0.0) {
                        return p2;
                    } else if (z < 0.0) {
                        return p1;
                    } else {
                        return math::Segment(p1, p2).project(X);
                    }
                } else if (y < 0.0) {
                    if (z < 0.0) {
                        return p0;
                    } else {
                        return math::Segment(p0, p2).project(X);
                    }
                } else if (z < 0.0) {
                    return math::Segment(p0, p1).project(X);
                }
                // we are in the triangle
                return X;
            }
            else {
                std::vector<Node> nodes = this->get<Node>();
                math::Point p0 = nodes[0].point();
                math::Point p1 = nodes[1].point();
                math::Point p2 = nodes[2].point();
                math::Point p3 = nodes[3].point();
                ts.emplace_back(p0,p1,p2);
                ts.emplace_back(p0,p2,p3);
                bool found_project=false;
                for(const auto& t:ts){
                    math::Point p0 = t.getPoint(0);
                    math::Point p1 = t.getPoint(1);
                    math::Point p2 = t.getPoint(2);
                    // we value the projected point
                    math::Point X = math::Plane(p0, p1, p2).project(AP);

                    TCoord x, y, z;
                    math::Point::computeBarycentric(p0, p1, p2, X, x, y, z);

                    if(x>=0 && y>=0 && z>=0) {
                        //we are in the triangle
                        found_project = true;


                        if (math::isZero(x, math::Constants::EPSILON)) {
                            if (math::isZero(y, math::Constants::EPSILON)) {
                                return p2;
                            } else if (math::isZero(z, math::Constants::EPSILON)) {
                                return p1;
                            } else {
                                return math::Segment(p1, p2).project(X);
                            }
                        } else if (math::isZero(y, math::Constants::EPSILON)) {
                            if (math::isZero(z, math::Constants::EPSILON)) {
                                return p0;
                            } else {
                                return math::Segment(p0, p2).project(X);
                            }
                        } else if (math::isZero(z, math::Constants::EPSILON)) {
                            return math::Segment(p0, p1).project(X);
                        }
                        // we are in the triangle
                        return X;
                    }
                }
                if(!found_project){
                    //means we are out of the square, we look for the closest projection on the quad
                    //boundary
                    std::vector<math::Point> candidates;
                    candidates.push_back(p1);
                    candidates.push_back(p2);
                    candidates.push_back(p3);
                    candidates.push_back(math::Segment(p0, p1).project(AP));
                    candidates.push_back(math::Segment(p1, p2).project(AP));
                    candidates.push_back(math::Segment(p2, p3).project(AP));
                    candidates.push_back(math::Segment(p3, p0).project(AP));
                    math::Point res = p0;
                    TCoord  dist = p0.distance2(AP);
                    for(const auto& c:candidates){
                        TCoord  dc = c.distance2(AP);
                        if(dc<dist){
                            dc = dist;
                            res = c;
                        }
                    }
                    return res;
                }

            }


    }
/*----------------------------------------------------------------------------*/
    void
    Face::getAdjacentNodes(const Node& ANode1, Node& ANode2, Node& ANode3)
    {
        TCellID id2, id3;
        getAdjacentNodes(ANode1.id(),id2,id3);
        ANode2 = m_owner->get<Node>(id2);
        ANode3 = m_owner->get<Node>(id3);
    }
    /*----------------------------------------------------------------------------*/
    void
    Face::getAdjacentNodes(const TCellID & ANode1, TCellID& ANode2, TCellID& ANode3)
    {
        std::vector<TCellID> nodes = getIDs<Node>();
        if (m_type == GMDS_QUAD) {
            int index1 = -1, index2, index3;
            for (int i = 0; i < 4; i++) {
                if (nodes[i] == ANode1)
                    index1 = i;
            }
            switch (index1) {
                case -1:
#ifdef _DEBUG
                    std::cout << "Index error" << std::endl;
#endif  //_DEBUG
                    throw GMDSException("getAdjacentNodes: node 1 is not adjacent to the face");
                    break;
                case 0:
                    index2 = 3;
                    index3 = 1;
                    break;
                case 1:
                    index2 = 0;
                    index3 = 2;
                    break;
                case 2:
                    index2 = 1;
                    index3 = 3;
                    break;
                case 3:
                    index2 = 2;
                    index3 = 0;
                    break;
            }

            ANode2 = nodes[index2];
            ANode3 = nodes[index3];
        } else if (m_type == GMDS_TRIANGLE) {
            int index1 = -1, index2, index3;
            for (int i = 0; i < 3; i++) {
                if (nodes[i] == ANode1)
                    index1 = i;
            }

            switch (index1) {
                case -1:
#ifdef _DEBUG
                    std::cout << "Index error" << std::endl;
                        for (int i = 0; i < 3; i++)
                                std::cout << nodes[i] << " ";
                        std::cout << "<- " << ANode1 << std::endl;
#endif  //_DEBUG
                    throw GMDSException("getAdjacentNodes: node 1 is not adjacent to the face");
                    break;
                case 0:
                    index2 = 2;
                    index3 = 1;
                    break;
                case 1:
                    index2 = 0;
                    index3 = 2;
                    break;
                case 2:
                    index2 = 1;
                    index3 = 0;
                    break;
                default:
                    throw GMDSException("Face::getAdjacentNodes index incorrectly computed.");
                    break;
            }

            ANode2 = nodes[index2];
            ANode3 = nodes[index3];
        } else if (m_type == GMDS_POLYGON) {
            bool found = false;
            unsigned int index1, index2, index3;
            for (unsigned int i = 0; !found; i++)
                if (nodes[i] == ANode1) {
                    found = true;
                    index1 = i;
                }

            const unsigned int last_index = nodes.size() - 1;

            if (!found) {
                std::cout << "Index error" << std::endl;
                throw GMDSException("getNodeOfOppositeEdge: node 1 is not adjacent to the face");
            } else if (index1 == 0) {
                index2 = last_index;
                index3 = 1;
            } else if (index1 == last_index) {
                index2 = last_index - 1;
                index3 = 0;
            } else {
                index2 = index1 - 1;
                index3 = index1 + 1;
            }

            ANode2 = nodes[index2];
            ANode3 = nodes[index3];
        } else
            throw GMDSException("getAdjacent node not availabe");
    }
/*----------------------------------------------------------------------------*/
    void
    Face::getOrderedEdges(std::vector<Edge>& AEdges) const
    {
        AEdges.clear();

            std::map<VirtualEdge, Edge> fakeEdgeMap;

            std::vector<Edge> edges = this->get<Edge>();
            for (auto & edge : edges) {
                    std::vector<Node> nodes = edge.get<Node>();
                    fakeEdgeMap[VirtualEdge(nodes[0].id(), nodes[1].id())] = edge;
            }

            std::vector<Node> nodes = this->get<Node>();
            for (unsigned int iNode = 0; iNode < nodes.size(); iNode++) {
                    if (fakeEdgeMap.find(VirtualEdge(nodes[iNode].id(), nodes[(iNode + 1) % nodes.size()].id())) ==
                        fakeEdgeMap.end()) {
                            throw GMDSException("Face::getOrderedEdges could not find edge");
                    }
                    AEdges.push_back(
                            fakeEdgeMap[VirtualEdge(nodes[iNode].id(), nodes[(iNode + 1) % nodes.size()].id())]);
            }
    }
/*----------------------------------------------------------------------------*/
    TCoord
    Face::distance(const math::Point& AP) const
    {
            return project(AP).distance(AP);
    }
/*----------------------------------------------------------------------------*/
    Face::Face(Mesh* AMesh, const ECellType AType, const TCellID& AID)
            : Cell(AMesh, AType, AID)  // mesh, type and id are filled in
    {
            //============================================
            // we keep a reference on the face container
            if (AMesh != nullptr) {
                    m_faces_container = AMesh->m_faces_container;
                    m_type_id = m_faces_container->getTypeID(AID);
            } else {
                    m_faces_container = nullptr;
                    m_type_id = NullID;
            }
    }
/*----------------------------------------------------------------------------*/
    Face::Face(const Face& AF)
            : Cell(AF.m_owner, AF.m_type, AF.m_id)
    {
            if (m_owner != nullptr)
                    m_faces_container = m_owner->m_faces_container;
            else
                    m_faces_container = nullptr;

            m_type_id = AF.m_type_id;
    }
/*----------------------------------------------------------------------------*/
    void
    Face::operator=(const Face& AF)
    {
            m_owner = AF.m_owner;
            m_type = AF.m_type;
            m_id = AF.m_id;
            if (m_owner != nullptr)
                    m_faces_container = m_owner->m_faces_container;
            else
                    m_faces_container = nullptr;
            m_type_id = AF.m_type_id;
    }
/*----------------------------------------------------------------------------*/
    bool
    Face::operator==(const Face& AFace) const
    {
            return (m_owner == AFace.m_owner && m_id == AFace.m_id);
    }
/*----------------------------------------------------------------------------*/
    bool
    Face::operator!=(const Face& AFace) const
    {
            return (!(*this == AFace));
    }
/*----------------------------------------------------------------------------*/
    Face::~Face()
    = default;
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGet(std::vector<Node>& ACells) const
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2N)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2N)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id].values(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_nodes_container->buildNode(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGet(std::vector<Edge>& ACells) const
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].values(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_edges_container->buildEdge(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGet(std::vector<Face>& ACells) const
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].values(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_faces_container->buildFace(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGet(std::vector<Region>& ACells) const
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].values(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].values(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_regions_container->buildRegion(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetNodeIDs(std::vector<TCellID>& ACells) const
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2N)[m_type_id].values(ACells);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2N)[m_type_id].values(ACells);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id].values(ACells);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetEdgeIDs(std::vector<TCellID>& ACells) const
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].values(ACells);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].values(ACells);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].values(ACells);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetFaceIDs(std::vector<TCellID>& ACells) const
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].values(ACells);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].values(ACells);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].values(ACells);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetRegionIDs(std::vector<TCellID>& ACells) const
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].values(ACells);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].values(ACells);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].values(ACells);
            }
    }

/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetAll(std::vector<Node>& ACells) const
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2N)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2N)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id].allValues(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_nodes_container->buildNode(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetAll(std::vector<Edge>& ACells) const
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].allValues(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_edges_container->buildEdge(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetAll(std::vector<Face>& ACells) const
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].allValues(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_faces_container->buildFace(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateGetAll(std::vector<Region>& ACells) const
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            std::vector<TCellID> cellIDs;
            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].allValues(cellIDs);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].allValues(cellIDs);
            } else
                    throw GMDSException("Not yet implemented");

            ACells.resize(cellIDs.size());
            for (unsigned int i = 0; i < cellIDs.size(); i++) {
                    ACells[i] = m_owner->m_regions_container->buildRegion(cellIDs[i]);
            }
    }
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2N))
                throw GMDSException("F2N adjacency is not supported by the mesh model");

        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2N)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2E))
                throw GMDSException("F2E adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2E)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2E)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2E)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2F))
                throw GMDSException("F2F adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2F)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2F)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2F)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
void Face::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
{
        if(!m_owner->m_model.has(F2R))
                throw GMDSException("F2R adjacency is not supported by the mesh model");
                
        if(m_type==GMDS_TRIANGLE){
                (*m_faces_container->m_T2R)[m_type_id].allValues(ACells);
        }       
        else if (m_type==GMDS_QUAD){
                (*m_faces_container->m_Q2R)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYGON){
                (*m_faces_container->m_P2R)[m_type_id].allValues(ACells);
        }
        else
                throw GMDSException("Not yet implemented");
}
/*----------------------------------------------------------------------------*/
    void
    Face::delegateSetNodeIDs(const std::vector<TCellID>& ACells)
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

        unsigned int nb_cells = nbNodes();
            if (m_type == GMDS_TRIANGLE) {
                    if (nb_cells != ACells.size()) {
                            std::cout << "OLD CELL NB: " << nb_cells << std::endl;
                            std::cout << "NEW CELL NB: " << ACells.size() << std::endl;
                            throw GMDSException("Invalid number of adj. entities");
                    }
                    (*m_faces_container->m_T2N)[m_type_id] = ACells;
            } else if (m_type == GMDS_QUAD) {
                    if (nb_cells != ACells.size())
                            throw GMDSException("Invalid number of adj. entities");

                    (*m_faces_container->m_Q2N)[m_type_id] = ACells;
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id] = ACells;
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateSetEdgeIDs(const std::vector<TCellID>& ACells)
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            unsigned int nb_cells = nbEdges();
            if (m_type == GMDS_TRIANGLE) {
                    if (nb_cells != ACells.size()) {
                        std::string mess =
                                "F2E for Triangle/ Invalid number of adj. entities: "+std::to_string(ACells.size()) +" instead of "+ std::to_string(nb_cells);
                        throw GMDSException(mess);
                    }

                    (*m_faces_container->m_T2E)[m_type_id] = ACells;
            } else if (m_type == GMDS_QUAD) {
                    if (nb_cells != ACells.size()){
                        std::string mess =
                                "F2E for Quad/ Invalid number of adj. entities: "+std::to_string(ACells.size()) +" instead of "+ std::to_string(nb_cells);
                        throw GMDSException(mess);
                    }

                (*m_faces_container->m_Q2E)[m_type_id] = ACells;
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id] = ACells;
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateSetFaceIDs(const std::vector<TCellID>& ACells)
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            unsigned int nb_cells = nbFaces();
            if (m_type == GMDS_TRIANGLE) {
                    if (nb_cells != ACells.size())
                            throw GMDSException("Invalid number of adj. entities");

                    (*m_faces_container->m_T2F)[m_type_id] = ACells;
            } else if (m_type == GMDS_QUAD) {
                    if (nb_cells != ACells.size())
                            throw GMDSException("Invalid number of adj. entities");

                    (*m_faces_container->m_Q2F)[m_type_id] = ACells;
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id] = ACells;
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
    {
        if (!m_owner->m_model.has(F2R))
            throw GMDSException("F2R adjacency is not supported by the mesh model");

        unsigned int nb_cells = nbRegions();
        if (nb_cells != ACells.size())
            throw GMDSException("Invalid number of adj. entities");

        if (m_type == GMDS_TRIANGLE) {
            (*m_faces_container->m_T2R)[m_type_id] = ACells;
        } else if (m_type == GMDS_QUAD) {
            (*m_faces_container->m_Q2R)[m_type_id] = ACells;
        } else if (m_type == GMDS_POLYGON) {
            (*m_faces_container->m_P2R)[m_type_id] = ACells;
        }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateNodeAdd(TCellID AID)
    {
        if (!m_owner->m_model.has(F2N))
            throw GMDSException("F2N adjacency is not supported by the mesh model");

        if (m_type == GMDS_TRIANGLE) {
            (*m_faces_container->m_T2N)[m_type_id].add(AID);
        } else if (m_type == GMDS_QUAD) {
            (*m_faces_container->m_Q2N)[m_type_id].add(AID);
        } else if (m_type == GMDS_POLYGON) {
            (*m_faces_container->m_P2N)[m_type_id].add(AID);
        }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateEdgeAdd(TCellID AID)
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].add(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].add(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].add(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateFaceAdd(TCellID AID)
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].add(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].add(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].add(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateRegionAdd(TCellID AID)
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].add(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].add(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].add(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateNodeRemove(TCellID AID)
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2N)[m_type_id].del(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2N)[m_type_id].del(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id].del(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateEdgeRemove(TCellID AID)
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].del(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].del(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].del(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateFaceRemove(TCellID AID)
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2F adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].del(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].del(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].del(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateRegionRemove(TCellID AID)
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].del(AID);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].del(AID);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].del(AID);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateNodeReplace(TCellID AID1, TCellID AID2)
    {
            if (!m_owner->m_model.has(F2N))
                    throw GMDSException("F2N adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2N)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2N)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2N)[m_type_id].replace(AID1, AID2);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateEdgeReplace(TCellID AID1, TCellID AID2)
    {
            if (!m_owner->m_model.has(F2E))
                    throw GMDSException("F2E adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2E)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2E)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2E)[m_type_id].replace(AID1, AID2);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateFaceReplace(TCellID AID1, TCellID AID2)
    {
            if (!m_owner->m_model.has(F2F))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2F)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2F)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2F)[m_type_id].replace(AID1, AID2);
            }
    }
/*----------------------------------------------------------------------------*/
    void
    Face::delegateRegionReplace(TCellID AID1, TCellID AID2)
    {
            if (!m_owner->m_model.has(F2R))
                    throw GMDSException("F2R adjacency is not supported by the mesh model");

            if (m_type == GMDS_TRIANGLE) {
                    (*m_faces_container->m_T2R)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_QUAD) {
                    (*m_faces_container->m_Q2R)[m_type_id].replace(AID1, AID2);
            } else if (m_type == GMDS_POLYGON) {
                    (*m_faces_container->m_P2R)[m_type_id].replace(AID1, AID2);
            }
    }
/*----------------------------------------------------------------------------*/
    std::ostream&
    operator<<(std::ostream& AStream, const Face& AF)
    {
            AStream << "Face " << AF.id() << " - ";
            if (AF.type() == GMDS_QUAD)
                    AStream << "Quad";
            else if (AF.type() == GMDS_TRIANGLE)
                    AStream << "Triangle";
            else if (AF.type() == GMDS_POLYGON)
                    AStream << "Polygon";
            else
                    AStream << "Invalid type";

            AStream << std::endl << "Nb Nodes " << AF.nbNodes();
            AStream << std::endl << "Nb Faces " << AF.nbFaces();
            return AStream;
    }
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
