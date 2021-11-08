/*----------------------------------------------------------------------------*/
#include <gmds/geodHoneyComb/RegularIcosahedron.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/utils/Array.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <map>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
RegularIcosahedron::
RegularIcosahedron(const math::Point &ACenter, const TCoord ARadius,
                   const int ANbSubdivision)
        : m_center(ACenter),
          m_radius(ARadius),
          m_subdivision(ANbSubdivision)
{
    m_representation=std::unique_ptr<Mesh>(new Mesh(MeshModel(DIM3|N|F|F2N|N2F)));

    // Nodes are simply located using the cartesian coordinates given in
    //  https://en.wikipedia.org/wiki/Regular_icosahedron
    TCoord psi = (1.0+ std::sqrt(5.0))/2.0;

    Node n1 = m_representation->newNode(0, 1, psi);
    Node n2 = m_representation->newNode(0, 1,-psi);
    Node n3 = m_representation->newNode(0,-1, psi);
    Node n4 = m_representation->newNode(0,-1,-psi);

    Node n5 = m_representation->newNode( 1, psi,0);
    Node n6 = m_representation->newNode( 1,-psi,0);
    Node n7 = m_representation->newNode(-1, psi,0);
    Node n8 = m_representation->newNode(-1,-psi,0);

    Node n9 = m_representation->newNode( psi,0, 1);
    Node n10= m_representation->newNode(-psi,0, 1);
    Node n11= m_representation->newNode( psi,0,-1);
    Node n12= m_representation->newNode(-psi,0,-1);

    m_representation->newTriangle(n1,n10,n3 );
    m_representation->newTriangle(n1,n3 ,n9 );
    m_representation->newTriangle(n1,n9 ,n5 );
    m_representation->newTriangle(n1,n5 ,n7 );
    m_representation->newTriangle(n1,n7 ,n10);

    m_representation->newTriangle(n5,n2,n7);
    m_representation->newTriangle(n7,n2,n12);
    m_representation->newTriangle(n7,n12,n10);
    m_representation->newTriangle(n10,n12,n8);
    m_representation->newTriangle(n10,n8,n3);
    m_representation->newTriangle(n3,n8,n6);
    m_representation->newTriangle(n3,n6,n9);
    m_representation->newTriangle(n9,n6,n11);
    m_representation->newTriangle(n9,n11,n5);
    m_representation->newTriangle(n5,n11,n2);

    m_representation->newTriangle(n11,n4,n2);
    m_representation->newTriangle(n2,n4,n12);
    m_representation->newTriangle(n12,n4,n8);
    m_representation->newTriangle(n8,n4,n6);
    m_representation->newTriangle(n6,n4,n11);

    for(auto n_id:m_representation->nodes()) {
        Node ni = m_representation->get<Node>(n_id);
        //first we translate the unit sphere to its center
        math::Point pi = ni.point();
        ni.setPoint(pi+m_center);
        //then we perform homothety
        pi = ni.point();
        TCoord  original_distance = pi.distance(m_center);
        TCoord ratio = m_radius/original_distance;
        math::Point pj = ratio*pi+(1-ratio)*m_center;
        ni.setPoint(pj);
    }

    MeshDoctor doc(reinterpret_cast<Mesh *>(&m_representation));
    doc.updateUpwardConnectivity();

    if(m_subdivision>1){
        subdivide();
    }
}
/*----------------------------------------------------------------------------*/
RegularIcosahedron::~RegularIcosahedron() {
}
/*----------------------------------------------------------------------------*/
void RegularIcosahedron::subdivide() {
    //we compute the location of new nodes on the edges.

    // new nodes per edge. They are computed from the edge end point with the
    // lowest id
    std::map<std::pair<TCellID ,TCellID >, std::vector<TCellID> > e2dn;
    for(auto f_id:m_representation->faces()){
        Face f = m_representation->get<Face>(f_id);
        std::vector<Node> f_nodes = f.get<Node>();
        for(auto i=0; i<3;i++){
            Node ni = f_nodes[i];
            Node nj = f_nodes[(i+1)%3];
            Node n_from = (ni.id()<nj.id())?ni:nj;
            Node n_to   = (ni.id()>nj.id())?ni:nj;
            if(e2dn.find(std::make_pair(n_from.id(),n_to.id()))==e2dn.end()){
                //not yet inserted
                std::vector<TCellID> edge_nodes;

                for(auto i=1; i<m_subdivision;i++){
                    TCoord w_from = (double)(m_subdivision-i)/(double)m_subdivision;
                    TCoord w_to = 1-w_from;
                    math::Point pij = w_from*n_from.point()+w_to*n_to.point();
                    pij = projectOnSphere(pij);
                    Node nij = m_representation->newNode(pij);
                    edge_nodes.push_back(nij.id());
                }
                e2dn[std::make_pair(n_from.id(),n_to.id())]= edge_nodes;
            }
        }
    }

    //now we build inner point in each face
    for(auto f_id:m_representation->faces()){
        Face f = m_representation->get<Face>(f_id);
        std::vector<Node> f_nodes = f.get<Node>();
        //starting from the first node of f, we get point of the boundary edges
        Node n0 = f_nodes[0];
        Node n1 = f_nodes[1];
        Node n2 = f_nodes[2];
        //We use an full 2D rectangular array while only a triangular one should
        // be better to use.
        TriArray<math::Point> grid(m_subdivision+1);
        TriArray<TCellID> ids(m_subdivision+1);
        //==================================================================
        //nO to n1 on the first line
        Node n_from = (n0.id()<n1.id())?n0:n1;
        Node n_to   = (n0.id()>n1.id())?n0:n1;
        bool invert =(n0.id()>n1.id());
        std::vector<TCellID> subdiv = e2dn[std::make_pair(n_from.id(), n_to.id())];
        if(invert)
            std::reverse(subdiv.begin(),subdiv.end());

        grid(0,0) = n0.point();
        ids(0,0) = n0.id();
        grid(m_subdivision,0) = n1.point();
        ids(m_subdivision,0) = n1.id();
        for(auto i=1; i<m_subdivision;i++){
            grid(i,0)= m_representation->get<Node>(subdiv[i-1]).point();
            ids(i,0) = m_representation->get<Node>(subdiv[i-1]).id();
        }
        //==================================================================
        //nO to n2 on the first column
        n_from = (n0.id()<n2.id())?n0:n2;
        n_to   = (n0.id()>n2.id())?n0:n2;
        invert =(n0.id()>n2.id());
        subdiv = e2dn[std::make_pair(n_from.id(), n_to.id())];
        if(invert)
            std::reverse(subdiv.begin(),subdiv.end());

        grid(0,m_subdivision) = n2.point();
        ids(0,m_subdivision) = n2.id();
        for(auto i=1; i<m_subdivision;i++){
            grid(0,i)= m_representation->get<Node>(subdiv[i-1]).point();
            ids(0,i) = m_representation->get<Node>(subdiv[i-1]).id();
        }
        //==================================================================
        //n2 to n1 on the first column
        n_from = (n2.id()<n1.id())?n2:n1;
        n_to   = (n2.id()>n1.id())?n2:n1;
        invert =(n2.id()>n1.id());
        subdiv = e2dn[std::make_pair(n_from.id(), n_to.id())];
        if(invert)
            std::reverse(subdiv.begin(),subdiv.end());

        for(auto i=1; i<m_subdivision;i++){
            grid(i,m_subdivision-i)= m_representation->get<Node>(subdiv[i-1]).point();
            ids(i,m_subdivision-i) = m_representation->get<Node>(subdiv[i-1]).id();
        }
        math::TransfiniteInterpolation::computeTri(grid);
        for(auto i=1;i<m_subdivision;i++){
            for(auto j=1;j<m_subdivision-i;j++){
                auto k = m_subdivision-i-j;
                math::Point pijk = projectOnSphere(grid(i,j,k));
                Node nijk = m_representation->newNode(pijk);
                ids(i,j,k) = nijk.id();
            }
        }
        m_representation->deleteFace(f_id);
        for(auto i=0;i<m_subdivision;i++){
            for(auto j=0;j<m_subdivision-i;j++){
                auto k = m_subdivision-i-j;
                m_representation->newTriangle(ids(i,j,k),
                                              ids(i+1,j,k-1),
                                              ids(i,j+1,k-1));
                if(k>1)
                    m_representation->newTriangle(ids(i+1,j+1,k-2),
                                                  ids(i+1,j,k-1),
                                                  ids(i,j+1,k-1));
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
math::Point RegularIcosahedron::projectOnSphere(const math::Point &APoint) {
    TCoord  original_distance = APoint.distance(m_center);
    TCoord ratio = m_radius/original_distance;
    return ratio*APoint+(1-ratio)*m_center;
}
/*----------------------------------------------------------------------------*/
std::unique_ptr<Mesh> RegularIcosahedron::getRepresentation()  {
    return std::move(m_representation);
}
/*----------------------------------------------------------------------------*/
void RegularIcosahedron::performQuadDualization() {
    std::map<TCellID , TCellID > f2dn;
    //create dual nodes of each face
    for(auto f_id:m_representation->faces()){
        Face f = m_representation->get<Face>(f_id);
        Node n = m_representation->newNode(projectOnSphere(f.center()));
        f2dn[f_id]=n.id();
    }

    //one node per edge
    std::map<std::pair<TCellID ,TCellID >, TCellID > e2dn;
    for(auto f_id:m_representation->faces()){
        Face f = m_representation->get<Face>(f_id);
        std::vector<Node> f_nodes = f.get<Node>();
        for(auto i=0; i<3;i++){
            Node ni = f_nodes[i];
            Node nj = f_nodes[(i+1)%3];
            if(e2dn.find(std::make_pair(ni.id(),nj.id()))==e2dn.end()){
                //not yet inserted
                math::Point pij = projectOnSphere(0.5*(ni.point()+nj.point()));
                Node nij = m_representation->newNode(pij);
                e2dn[std::make_pair(ni.id(),nj.id())]=nij.id();
            }
        }
    }


    //for each vertex of the icosahedron (index 0 to 11), we store the id of adjacent
    // vertices and faces in a circular order
    TCellID tab_dual_node[12][5]={
            {2,9,6,4,8}, //vertex 0
            {6,11,3,10,4}, //vertex 1
            {9,0,8,5,7}, //vertex 2
            {1,11,7,5,10}, //vertex 3
            {8,0,6,1,10}, //vertex 4
            {2,8,10,3,7}, //vertex 5
            {0,9,11,1,4}, //vertex 6
            {2,5,3,11,9}, //vertex 7
            {0,4,10,5,2}, //vertex 8
            {0,2,7,11,6},  //vertex 9
            {1,3,5,8,4},  //vertex 10
            {1,6,9,7,3} //vertex 11
    };

    TCellID tab_dual_face[12][5]={
            {0,4,3,2,1}, //vertex 0
            {6,16,15,14,5}, //vertex 1
            {0,1,11,10,9}, //vertex 2
            {16,17,18,19,15}, //vertex 3
            {2,3,5,14,13}, //vertex 4
            {11,12,19,18,10}, //vertex 5
            {4,7,6,5,3}, //vertex 6
            {10,18,17,8,9}, //vertex 7
            {2,13,12,11,1}, //vertex 8
            {0,9,8,7,4},  //vertex 9
            {15,19,12,13,14},  //vertex 10
            {6,7,8,17,16}//vertex 11
    };

    for(auto f_id:m_representation->faces()){
        m_representation->deleteFace(f_id);
    }
    for(auto i=0;i<12;i++){
        //we create five quads
        TCellID i0 = i; //first node id
        for(auto j=0;j<5;j++){
            TCellID i1 = tab_dual_node[i][j];
            TCellID i3 = tab_dual_node[i][(j+1)%5];
            TCellID i2 = tab_dual_face[i][j];
            m_representation->newQuad(i0,
                                      e2dn[std::make_pair(i0,i1)],
                                      f2dn[i2],
                                      e2dn[std::make_pair(i0,i3)]);
        }
    }
}
/*----------------------------------------------------------------------------*/
