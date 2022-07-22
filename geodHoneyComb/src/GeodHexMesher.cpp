/*----------------------------------------------------------------------------*/
#include <gmds/geodHoneyComb/GeodHexMesher.h>
#include <gmds/geodHoneyComb/RegularIcosahedron.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
GeodHexMesher::GeodHexMesher()
: m_radius(0), m_center(0,0,0),m_is_sphere(GEOD_NOT_DEFINED)
{}
/*----------------------------------------------------------------------------*/
GeodHexMesher::~GeodHexMesher() {}
/*----------------------------------------------------------------------------*/
void GeodHexMesher::setRadius(const double AR) {
    m_radius=AR;
}
/*----------------------------------------------------------------------------*/
void GeodHexMesher::setCenter(const math::Point &AP) {
    m_center=AP;
}
/*----------------------------------------------------------------------------*/
void GeodHexMesher::
setLayerData(std::map<double, math::DiscretizationScheme1D *> &ALayerData) {
    for (auto data: ALayerData) {
        if (data.first < 0)
            throw GMDSException("Negative radius value for a layer");
        if (data.first >= m_radius)
            throw GMDSException("Maximum layer radius cannot be greater or equal to the sphere radius");
    }
    m_layer_data=ALayerData;
}

/*----------------------------------------------------------------------------*/
std::unique_ptr<Mesh> GeodHexMesher::getMesh()  {
    return std::move(m_mesh);
}
/*----------------------------------------------------------------------------*/
GeodHexMesher::OpResult GeodHexMesher::execute() {
    // we first check that the data provided are coherent (forgotten setter or
    // several set that could have led to inconsistency
    if(m_radius<=0){
        return GEOD_FAILURE_WRONG_RADIUS;
    }
    m_is_sphere=false;
    for (auto data: m_layer_data) {
        if(data.first==0) //we found a layer starting from the geode center
            m_is_sphere=true;
        if (data.first < 0)
            return GEOD_FAILURE_LAYER_RADIUS_NEGATIVE;
        if (data.first >= m_radius)
            return GEOD_FAILURE_LAYER_RADIUS_GREATER_THAN_SPHERE_RADIUS;
    }

    //Now we can create the mesh!!!!
    RegularIcosahedron ico(m_center,m_radius,5);
    ico.performQuadDualization();

     m_mesh= ico.getRepresentation();
    //surface mesh has the following model DIM3|N|F|F2N|N2F
    //we update the model to add R and R2N
    MeshModel mod = m_mesh->getModel();
    mod.add(R|R2N);
    m_mesh->changeModel(mod);

    //===================================================================
    // we mark all the surface nodes
    TInt sphere_mark = m_mesh->newMark<Node>();
    for(auto n_id:m_mesh->nodes()){
        m_mesh->mark<Node>(n_id,sphere_mark);
    }

    //===================================================================
    // We create inner nodes, and we store for each existing sphere node
    // the line of nodes coming from the surface to the center
    TCellID sphere_center_node_id=NullID;
    if(m_is_sphere){
        Node n_center=m_mesh->newNode(m_center);
        sphere_center_node_id=n_center.id();
    }
    std::map<TCellID ,std::vector<TCellID> > surf_to_center_line;
    std::map<TCellID ,math::DiscretizationScheme1D*> surf_to_center_discretization;

    for(auto layer:m_layer_data){
        double layer_radius = layer.first;
        for(auto n_id:m_mesh->nodes()){
            Node ni = m_mesh->get<Node>(n_id);
            if(m_mesh->isMarked(ni,sphere_mark)){
                //means the node is a surface one and not a new line node
                if(layer_radius==0){
                    //center point which is already created and shared
                    //by all the lines
                    surf_to_center_line[n_id].push_back(sphere_center_node_id);
                    surf_to_center_discretization[sphere_center_node_id] = layer.second;
                }
                else{
                    //we create the node
                    math::Vector v=ni.point()-m_center;
                    v.normalize();
                    v = layer_radius*v;
                    Node nj = m_mesh->newNode(m_center+v);
                    surf_to_center_line[n_id].push_back(nj.id());
                    surf_to_center_discretization[nj.id()] = layer.second;
                }
            }
        }
    }
    //we initialize the line with their head, i.e. the surface node that already exists
    for (auto n_id: m_mesh->nodes()) {
        if (m_mesh->isMarked<Node>(n_id, sphere_mark)) {
            surf_to_center_line[n_id].push_back(n_id);
        }
    }
    // lines were created from center to surface, we invert them so
    for(auto &line:surf_to_center_line){
        std::reverse(line.second.begin(),line.second.end());
    }
    //===================================================================
    // We create inner layer nodes and we store them in a map where the
    // input is the id of the node with the smallest radius.
    std::map<TCellID ,std::vector<TCellID> > inner_layer_nodes;
  //  bool first_line=true;
    for(auto line_data:surf_to_center_line){
        std::vector<TCellID> line_nodes = line_data.second;

        // line_nodes contains the ordered list of block nodes for the
        // line going from the surface to the center node.
        for(auto i=0;i<line_nodes.size()-1;i++){
            Node origin_node = m_mesh->get<Node>(line_nodes[i]);
            math::Point orig_pnt =origin_node.point();
            Node dest_node = m_mesh->get<Node>(line_nodes[i+1]);
            math::Point dest_pnt = dest_node.point();
  //         if (first_line)
    //           std::cout<<"from "<<orig_pnt<<" to "<<dest_pnt<<std::endl;

            math::DiscretizationScheme1D* d = surf_to_center_discretization[dest_node.id()];
            d->setOrigin(orig_pnt);
            d->setDestination(dest_pnt);
            int nb_pnts = d->getNbPoints();
            //we store the origin node id to avoid conditional afterward
            inner_layer_nodes[line_nodes[i]].push_back(origin_node.id());

            for(auto j=1;j<nb_pnts-1;j++){
                Node nj = m_mesh->newNode((*d)(j));
//                if (first_line)
//                    std::cout<<"add "<<j<<" -> "<<nj.point()<<std::endl;
                inner_layer_nodes[line_nodes[i]].push_back(nj.id());
            }
            //we store the destination node id to avoid conditional afterward
            inner_layer_nodes[line_nodes[i]].push_back(dest_node.id());
        }
    //    first_line=false;

    }
    //===================================================================
    // all nodes are created, we can now build 3D cells
    // we create one variable with the layer for regions
    Variable<int>* layer_var = m_mesh->newVariable<int,GMDS_REGION>("layer_index");
    // now for each face of the mesh (i.e. on the sphere) we create a line of region
    // cells
    for(auto f_id:m_mesh->faces()){
        Face f = m_mesh->get<Face>(f_id);
        std::vector<Node> f_nodes = f.get<Node>();
        if(f_nodes.size()!=4){
            throw GMDSException("Geode surface must be made of quad faces only!");
        }
        std::vector<TCellID> line0 = surf_to_center_line[f_nodes[0].id()];
        std::vector<TCellID> line1 = surf_to_center_line[f_nodes[1].id()];
        std::vector<TCellID> line2 = surf_to_center_line[f_nodes[2].id()];
        std::vector<TCellID> line3 = surf_to_center_line[f_nodes[3].id()];

        for(auto i=0;i<line0.size()-1;i++){
            //we are in a layer now
            std::vector<TCellID> l0 = inner_layer_nodes[line0[i]];
            std::vector<TCellID> l1 = inner_layer_nodes[line1[i]];
            std::vector<TCellID> l2 = inner_layer_nodes[line2[i]];
            std::vector<TCellID> l3 = inner_layer_nodes[line3[i]];
            for(auto j=0;j<l0.size()-1;j++){
                Region r;
                if(i==line0.size()-2 && j==l0.size()-2 && m_is_sphere){
                    r= m_mesh->newPyramid(l0[j],
                                          l1[j],
                                          l2[j],
                                          l3[j],
                                          l0[j+1]);
                }
                else{
                    r= m_mesh->newHex(l0[j+1],
                                      l1[j+1],
                                      l2[j+1],
                                      l3[j+1],
                                      l0[j],
                                      l1[j],
                                      l2[j],
                                      l3[j]);

                }
                layer_var->set(r.id(),i+1);
            }
        }
    }
    return GeodHexMesher::GEOD_SUCCESS;
}
/*----------------------------------------------------------------------------*/
