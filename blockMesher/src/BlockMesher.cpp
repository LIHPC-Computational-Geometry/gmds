/*----------------------------------------------------------------------------*/
#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/utils/LocalCellTopology.h>
#include <gmds/math/TransfiniteInterpolation.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BlockMesher::BlockMesher(Mesh *ABlocks, cad::GeomMeshLinker *ALinker) :
        m_blocks(ABlocks),
        m_linker(ALinker)
{
    m_with_geometry=(ALinker!=NULL);
    if(ALinker->mesh()!=ABlocks)
        throw GMDSException("Invalid linker provided to a BlockMesher instance");
    m_nb_edges_in_discretization=0;
    m_mesh = new Mesh(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N));
}
/*----------------------------------------------------------------------------*/
BlockMesher::~BlockMesher() {
    delete m_mesh;
}
/*----------------------------------------------------------------------------*/
BlockMesher::STATUS BlockMesher::execute(const int ANb) {
    m_nb_edges_in_discretization=ANb;
    //==========================================================================
    // 1. We mesh block vertices
    if(!meshVertices())
        return BlockMesher::BLOCK_VERTEX_MESHING_ERROR;
    //==========================================================================
    // 2. We mesh block edges

    if(!meshEdges())
        return BlockMesher::BLOCK_EDGE_MESHING_ERROR;
    //==========================================================================
    // 3. We mesh block faces

    if(!meshFaces())
        return BlockMesher::BLOCK_FACE_MESHING_ERROR;
    //==========================================================================
    // 4. We mesh blocks

    //if(!meshBlocks())
    //    return BlockMesher::BLOCK_MODEL_ERROR;

    return BlockMesher::SUCCESS;
}
/*----------------------------------------------------------------------------*/
BlockMesher::STATUS BlockMesher::execute(const double AEdgeTargetSize) {
    return BlockMesher::NOT_YET_IMPLEMENTED;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher::meshVertices() {
    for (auto v_id : m_blocks->nodes()){
        Node v = m_blocks->get<Node>(v_id);
        Node n = m_mesh->newNode(v.getPoint());
        m_blockV_to_meshN[v.id()]=n.id();
        if(m_with_geometry){
            Cell::Data d(0, m_linker->getGeomId(v));
            m_mesh_node_classification[n.id()]=d;
        }
    }
    return true;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher::meshEdges() {

    for (auto e_id : m_blocks->edges()) {
        Edge e = m_blocks->get<Edge>(e_id);
        auto n_ids = e.getIDs<Node>();
        if (n_ids.size() != 2)
            return false;

        //the edge must be discritize in the gmds id increasing direction
        TCellID id_node_0 = NullID;
        TCellID id_node_1 = NullID;
        if (n_ids[0] < n_ids[1]) {
            id_node_0 = n_ids[0];
            id_node_1 = n_ids[1];
        } else if (n_ids[0] > n_ids[1]) {
            id_node_0 = n_ids[1];
            id_node_1 = n_ids[0];
        } else {
            //the edge ends are the same
            return false;
        }
        Node mesh_n0 = m_mesh->get<Node>(m_blockV_to_meshN[id_node_0]);
        Node mesh_n1 = m_mesh->get<Node>(m_blockV_to_meshN[id_node_1]);
        math::Point p0 = mesh_n0.getPoint();
        math::Point p1 = mesh_n1.getPoint();
        auto N = m_nb_edges_in_discretization;
        std::vector<TCellID> discretization_ids;
        discretization_ids.resize(m_nb_edges_in_discretization + 1);
        discretization_ids[0]=mesh_n0.id();
        double inv = 1.0/(double)(N);
        for (auto i = 1; i < m_nb_edges_in_discretization; i++) {
            math::Point pi = inv * ((N - i) * p0 + i * p1);
            Node ni = m_mesh->newNode(pi);
            discretization_ids[i] = ni.id();
        }
        discretization_ids[m_nb_edges_in_discretization]=mesh_n1.id();
        //we store the edge discretization
        m_blockE_to_meshN.insert(std::make_pair(e.id(),discretization_ids));


        if(m_with_geometry){
            //all nodes internal to the block edges are classified as the edge
            for (auto i = 1; i < m_nb_edges_in_discretization; i++) {
                if(m_linker->getGeomDim(e)==cad::GeomMeshLinker::LINK_CURVE) {
                    m_mesh_node_classification[discretization_ids[i]] = Cell::Data(1, m_linker->getGeomId(e));
                } else  if(m_linker->getGeomDim(e)==cad::GeomMeshLinker::LINK_SURFACE) {
                    m_mesh_node_classification[discretization_ids[i]] = Cell::Data(2, m_linker->getGeomId(e));
                }
            }
            //now, we create the edge and link them
            if(m_linker->getGeomDim(e)==cad::GeomMeshLinker::LINK_CURVE){
                //we create mesh edges and we classified them on the curve!
                for (auto i = 1; i <= m_nb_edges_in_discretization; i++) {
                    Edge e_i = m_mesh->newEdge(discretization_ids[i-1],discretization_ids[i]);
                    Cell::Data d(1, m_linker->getGeomId(e));
                    m_mesh_edge_classification[e_i.id()]=d;
                }
            }
        }
    }
    return true;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher::
getEdgeFrom(const TCellID AN0, const TCellID AN1,
            const std::vector<Edge> &AEdges, Edge &AResult,
            bool &AInverted){
    for(auto e: AEdges){
        std::vector<TCellID > e_node_ids = e.getIDs<Node>();
        if(e_node_ids[0]==AN0 && e_node_ids[1]==AN1){
            AResult=e;
            AInverted=false;
            return true;
        } else if(e_node_ids[0]==AN1 && e_node_ids[1]==AN0){
            AResult=e;
            AInverted=true;
            return true;
        }
    }
    return false;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher::meshFaces()  {

    for(auto f_id:m_blocks->faces()){
        Face f = m_blocks->get<Face> (f_id);
        std::vector<Edge> f_edges = f.get<Edge>();
        std::vector<Node> f_nodes = f.get<Node>();
        //Edge locally defined from 0 to 1
        Node n0 = f_nodes[0];
        Node n1 = f_nodes[1];
        bool inverted_01 = false;
        Edge e01;
        if(!getEdgeFrom(n0.id(),n1.id(), f_edges, e01,inverted_01))
            return false;

        //Edge locally defined from 1 to 2
        Node n2 = f_nodes[2];
        bool inverted_12 = false;
        Edge e12;
        if(!getEdgeFrom(n1.id(),n2.id(), f_edges, e12,inverted_12))
            return false;

        //Edge locally defined from 3 to 2 (to be in the same direction as e01
        Node n3 = f_nodes[3];
        bool inverted_32 = false;
        Edge e32;
        if(!getEdgeFrom(n3.id(),n2.id(), f_edges, e32,inverted_32))
            return false;

        //Edge locally defined from 0 to 3 (to be in the same direction as e12
        bool inverted_03 = false;
        Edge e03;
        if(!getEdgeFrom(n0.id(),n3.id(), f_edges, e03,inverted_03))
            return false;

        //now we get the 4 edges in the expected direction. We can store
        //their respective discretization in the right way into a 2D-like array
        std::vector<std::vector<TCellID> > discretization_ids;
        std::vector<std::vector<math::Point> > discretization_pnts;
        discretization_ids.resize(m_nb_edges_in_discretization+1);
        for(auto &v:discretization_ids){
            v.resize(m_nb_edges_in_discretization+1);
        }
        discretization_pnts.resize(m_nb_edges_in_discretization+1);
        for(auto &v:discretization_pnts){
            v.resize(m_nb_edges_in_discretization+1);
        }
        auto N = m_nb_edges_in_discretization+1;
        // Corner points retrieval
        discretization_ids[0  ][0  ]= n0.id();
        discretization_ids[N-1][0  ]= n1.id();
        discretization_ids[N-1][N-1]= n2.id();
        discretization_ids[0  ][N-1]= n3.id();
        discretization_pnts[0  ][0  ]= n0.getPoint();
        discretization_pnts[N-1][0  ]= n1.getPoint();
        discretization_pnts[N-1][N-1]= n2.getPoint();
        discretization_pnts[0  ][N-1]= n3.getPoint();

        //Boundary edges point retrieval for edge 01
        std::vector<TCellID> edge_disc = m_blockE_to_meshN[e01.id()];

        if(n0.id()!=edge_disc[0])
           std::reverse(edge_disc.begin(), edge_disc.end());
        for(auto i=1; i<N-1;i++){
            discretization_ids[i][0]= edge_disc[i];
            discretization_pnts[i][0]=m_mesh->get<Node>(edge_disc[i]).getPoint();
        }
        //Boundary edges point retrieval for edge 01
        edge_disc = m_blockE_to_meshN[e12.id()];
        if(n1.id()!=edge_disc[0])
            std::reverse(edge_disc.begin(), edge_disc.end());

        for(auto i=1; i<N-1;i++){
            discretization_ids[N-1][i]= edge_disc[i];
            discretization_pnts[N-1][i]=m_mesh->get<Node>(edge_disc[i]).getPoint();
        }
        //Boundary edges point retrieval for edge 32
        edge_disc = m_blockE_to_meshN[e32.id()];
        if(n3.id()!=edge_disc[0])
            std::reverse(edge_disc.begin(), edge_disc.end());
        for(auto i=1; i<N-1;i++){
            discretization_ids[i][N-1]= edge_disc[i];
            discretization_pnts[i][N-1]=m_mesh->get<Node>(edge_disc[i]).getPoint();
        }
        //Boundary edges point retrieval for edge 01
        edge_disc = m_blockE_to_meshN[e03.id()];
        if(n0.id()!=edge_disc[0])
            std::reverse(edge_disc.begin(), edge_disc.end());
        for(auto i=1; i<N-1;i++){
            discretization_ids[0][i]= edge_disc[i];
            discretization_pnts[0][i]=m_mesh->get<Node>(edge_disc[i]).getPoint();
        }
        //compute the points location and create node
        math::TransfiniteInterpolation::compute(discretization_pnts);
        for(auto i=1; i<N-1;i++) {
            for(auto j=1; j<N-1;j++) {
                Node nij = m_mesh->newNode(discretization_pnts[i][j]);
                discretization_ids[i][j]= nij.id();
                if(m_with_geometry && m_linker->getGeomDim(f)==cad::GeomMeshLinker::LINK_SURFACE){
                    //classified on a surface and so the node too
                    m_mesh_node_classification[nij.id()]=Cell::Data(2, m_linker->getGeomId(f));

                }
            }
        }
        m_blockF_to_meshN[f.id()]=discretization_ids;
    }

    Variable<int>* var_mesh_surf = m_mesh->getOrCreateVariable<int,GMDS_FACE>("surf_color");
    for(auto f_id:m_blocks->faces()){
        std::vector<TCellID> face_ids;
        std::vector<std::vector<TCellID> > n_ids = m_blockF_to_meshN[f_id];
        auto I = n_ids.size();
        auto J = n_ids[0].size();
        for(auto i=0;i<I-1;i++){
            for(auto j=0;j<J-1;j++){
                Face q = m_mesh->newQuad(n_ids[i][j],n_ids[i+1][j],n_ids[i+1][j+1],n_ids[i][j+1]);
                face_ids.push_back(q.id());
                var_mesh_surf->value(q.id())=f_id;
                if(m_with_geometry && m_linker->getGeomDim<Face>(f_id)==cad::GeomMeshLinker::GeomMeshLinker::LINK_SURFACE){
                    //classified on a surface and so the node too
                    m_mesh_face_classification[q.id()]=Cell::Data(2, m_linker->getGeomId<Face>(f_id));
                }

            }
        }
    }
    return true;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher::meshBlocks()  {
    return false;
}

/*----------------------------------------------------------------------------*/
