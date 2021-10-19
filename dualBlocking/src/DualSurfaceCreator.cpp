//
// Created by calderans on 07/03/19.
//


#include <gmds/math/Numerics.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include "gmds/dualBlocking/DualSurfaceCreator.h"
#include "Predicates_psm.h"

using namespace gmds;
using namespace db;

/*----------------------------------------------------------------------------*/

DualSurfaceCreator::DualSurfaceCreator(Mesh *AMesh,
        //TCellID AID,
        //int ASheetID,
                                       gmds::Variable<std::vector<intersectInfo>>* info_Var,
                                       gmds::Variable<gmds::math::Chart>* vertex_chart_Var,
                                       gmds::Variable<int>* propagation_round_Var,
                                       gmds::Variable<int>* singularities_Var,
                                       int wave_precedent_mark,
                                       int face_treated_mark,
                                       int wave_tet_mark,
                                       double AMin_length,
                                       cad::GeomMeshLinker &ALinker,
                                       cad::FACManager &AManager):
        DualSheetCreator(AMesh,/* AID, ASheetID,*/ info_Var, vertex_chart_Var, propagation_round_Var, singularities_Var, wave_precedent_mark,face_treated_mark,wave_tet_mark,AMin_length,ALinker,AManager) {

    MeshDoctor doci(&m_smesh);
    doci.buildFacesAndR2F();
    doci.buildEdgesAndX2E();
    doci.updateUpwardConnectivity();
    wave_id = m_smesh.newVariable<int, GMDS_FACE>("wave");
    edge_Gdim = m_smesh.newVariable<int, GMDS_EDGE>("edge_dim");
    node_Gdim = m_smesh.newVariable<int, GMDS_NODE>("node_dim");

    try {
        m_part = m_mesh->getVariable<int,GMDS_REGION>("CUT_tet");
    }catch (GMDSException e){

        m_part = m_mesh->newVariable<int,GMDS_REGION>("CUT_tet");
        for(auto r : m_mesh->regions()){
            m_part->set(r,1);
        }
    }

}
/*----------------------------------------------------------------------------*/

DualSurfaceCreator::~DualSurfaceCreator(){

}
/*----------------------------------------------------------------------------*/

void DualSurfaceCreator::setNorm(math::Vector3d const &norm) {
    start_norm = norm;
}
/*----------------------------------------------------------------------------*/

gmds::math::Vector3d DualSurfaceCreator::getNorm() {
    return start_norm;
}
/*----------------------------------------------------------------------------*/

bool DualSurfaceCreator::execute() {

    propagation_round = 0;

    surface.clear();

    for (auto r : m_mesh->regions()) {
        (*m_propagation_round)[r] = -1;
    }


    //bool result = propagation_loop(tetID, start_norm);
    //return result;

    //============================================================
    //Initialisation step
    //============================================================

    Region current_region = m_mesh->get<Region>(tetID);
    Edge current_edge;
    math::Vector3d norm;

    math::Point plane_point;

    plane_point = current_region.center();

    norm = findNormal(plane_point, start_norm, current_region);
    std::vector<math::Vector3d> norm_vec;
    norm_vec.push_back(norm);
    surface.emplace(std::pair<TCellID, std::vector<math::Vector3d>>(tetID, norm_vec));

    (*m_propagation_round)[tetID] = 0;

    math::Plane plane(plane_point, norm);

    std::cout << "norm = " << norm.X() << "," << norm.Y() << "," << norm.Z() << std::endl;

    bool moved_node_init = false;
    TCellID moved_node_id_init = NullID;
    std::vector<intersectInfo> infos = intersect(tetID, plane, moved_node_init, moved_node_id_init);

    while (moved_node_init) {

        //If the plan in tetID is on node, we randomly move the node and repeat

        moved_node_init = false;

        Node node_moved = m_mesh->get<Node>(moved_node_id_init);

        std::vector<Edge> adj_edges = node_moved.get<Edge>();
        double d = 10000;
        for (auto const &a: adj_edges) {
            std::vector<Node> a_nodes = a.get<Node>();
            Node opp_node = (a_nodes[0].id() == moved_node_id_init) ? a_nodes[1] : a_nodes[0];
            double ad = opp_node.point().distance(node_moved.point());
            if (ad < d)
                d = ad;
        }

        d = d / (5 * sqrt(3));
        double dx = drand48() * d;
        double dy = drand48() * d;
        double dz = drand48() * d;

        math::Point new_loc(node_moved.point().X() + dx,
                            node_moved.point().Y() + dy,
                            node_moved.point().Z() + dz);
        node_moved.setPoint(new_loc);


        infos = intersect(tetID, plane, moved_node_init, moved_node_id_init);

    }
    moved_node_init = false;
    std::vector<math::Vector3d> infos_vectors;


    for (unsigned int i = 0; i < infos.size(); i++) {
        infos_vectors.push_back(infos[i].v);
    }


    //recomputing the normal plan and repeating the cut
    norm = meanVector(infos_vectors);
    plane.set(plane_point, norm);

    infos.clear();
    infos = intersect(tetID, plane, moved_node_init, moved_node_id_init);
    if (moved_node_init) {
        infos = intersect(tetID, plane, moved_node_init, moved_node_id_init);
    }

    //mark the edges of the first tet intersected as previous edges
    for (auto const &i : infos) {
        (*var_info)[i.e].push_back(i);
        m_mesh->mark(m_mesh->get<Edge>(i.e), wave_precedent);
    }


    propagation_round = 1;
    std::set<TCellID> wave_tet, next_tet;

    //here we create the first wave by getting tet from all intersect edges
    for (auto const &iI : infos) {
        std::vector<TCellID> tets = m_mesh->get<Edge>(iI.e).getIDs<Region>();
        for (auto t : tets) {
            if ((*m_propagation_round)[t] == -1) {
                if (!isAxisSingularity(t, iI.v) && m_part->value(t) <= 1) {
                    wave_tet.insert(t);
                    (*m_propagation_round)[t] = propagation_round;
                } else {
                    std::cout<<"Return 2?"<<std::endl;
                    return false;
                }
            }
        }
    }


    infos.clear();

    bool result = propagation_loop(wave_tet);


    std::cout<<"Surface size = "<<surface.size()<<std::endl;
    if(!result){
        //surface.clear();
    }

    //clear();

    return result;

}
/*----------------------------------------------------------------------------*/
gmds::Mesh* DualSurfaceCreator::buildSurfaceSheet() {
    // TEST: Fonction qui affiche la surface duale

    std::cout<<"Plan surface creation"<<std::endl;

    m_smesh.clear();

    std::map<int,Node> edge_to_node;

    for(auto e : m_mesh->edges()){
        if(!(*var_info)[e].empty()){
            Node n = m_smesh.newNode((*var_info)[e][0].p);
            edge_to_node.emplace(e,n);
            (*node_Gdim)[n.id()] = linker.getGeomDim<Edge>(e);
        }
    }

    for(auto t : surface){


        Region r = m_mesh->get<Region>(t.first);
        std::vector<Node> nodes;

        for(auto e : r.getIDs<Edge>()){
            if(!(*var_info)[e].empty()){
                nodes.push_back(edge_to_node[e]);
            }
        }

        if(nodes.size() > 2){
            if(nodes.size() == 4){
                math::Triangle tri1(nodes[0].point(),nodes[1].point(),nodes[2].point());
                math::Triangle tri2(nodes[0].point(),nodes[2].point(),nodes[3].point());

                math::Vector3d vec1 = tri1.getNormal();
                math::Vector3d vec2 = tri2.getNormal();

                if(vec1.dot(vec2)< 0){
                    Node tmp = nodes[3];
                    nodes[3] = nodes[2];
                    nodes[2] = tmp;
                }
                math::Triangle tri3(nodes[0].point(),nodes[1].point(),nodes[2].point());
                math::Triangle tri4(nodes[0].point(),nodes[2].point(),nodes[3].point());

                math::Vector3d vec3 = tri3.getNormal();
                math::Vector3d vec4 = tri4.getNormal();
                if(vec3.dot(vec4)< 0){
                    Node tmp = nodes[2];
                    nodes[2] = nodes[1];
                    nodes[1] = tmp;
                }
            }

            Face f = m_smesh.newFace(nodes);
            (*wave_id)[f.id()] = (*m_propagation_round)[t.first];

            Edge temoin;
            if(nodes.size() == 4){
                bool exist = false;
                for(auto e : nodes[0].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[0] && e.get<Node>()[1] == nodes[1]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[0] && e.get<Node>()[0] == nodes[1]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[0],nodes[1]);
                    nodes[0].add<Edge>(e_new);
                    nodes[1].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;

                for(auto e : nodes[1].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[1] && e.get<Node>()[1] == nodes[2]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[1] && e.get<Node>()[0] == nodes[2]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[1],nodes[2]);
                    nodes[1].add<Edge>(e_new);
                    nodes[2].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;

                for(auto e : nodes[2].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[2] && e.get<Node>()[1] == nodes[3]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[2] && e.get<Node>()[0] == nodes[3]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[2],nodes[3]);
                    nodes[2].add<Edge>(e_new);
                    nodes[3].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;

                for(auto e : nodes[3].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[3] && e.get<Node>()[1] == nodes[0]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[3] && e.get<Node>()[0] == nodes[0]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[3],nodes[0]);
                    nodes[3].add<Edge>(e_new);
                    nodes[0].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;
            }

            if(nodes.size() == 3){
                bool exist = false;
                for(auto e : nodes[0].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[0] && e.get<Node>()[1] == nodes[1]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[0] && e.get<Node>()[0] == nodes[1]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[0],nodes[1]);
                    nodes[0].add<Edge>(e_new);
                    nodes[1].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;

                for(auto e : nodes[1].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[1] && e.get<Node>()[1] == nodes[2]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[1] && e.get<Node>()[0] == nodes[2]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[1],nodes[2]);
                    nodes[1].add<Edge>(e_new);
                    nodes[2].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;


                for(auto e : nodes[2].get<Edge>()){
                    if(e.get<Node>()[0] == nodes[2] && e.get<Node>()[1] == nodes[0]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }else if(e.get<Node>()[1] == nodes[2] && e.get<Node>()[0] == nodes[0]){
                        e.add<Face>(f);
                        exist = true;
                        break;
                    }
                }
                if(!exist){
                    Edge e_new = m_smesh.newEdge(nodes[2],nodes[0]);
                    nodes[2].add<Edge>(e_new);
                    nodes[0].add<Edge>(e_new);
                    e_new.add<Face>(f);
                }
                exist = false;
            }
        }
    }
    //std::cout<<"Nb edges = "<<m_smesh.getNbEdges()<<std::endl;

    int cpt_out = 0;
    int cpt_in = 0;

    for(auto f : m_smesh.faces()){
        Face f_ = m_smesh.get<Face>(f);
        if(f_.get<Node>().size() == 2){
            std::cout<<"Face 2 nodes"<<std::endl;
            exit(2);
        }
    }


    MeshDoctor doci(&m_smesh);
    doci.buildE();
    doci.updateUpwardConnectivity();
    //std::cout<<"Nb edges = "<<m_smesh.getNbEdges()<<std::endl;
    //Loop edge classification
    for(auto eid : m_smesh.edges()){


        //std::cout<<"Test"<<std::endl;
        Edge e = m_smesh.get<Edge>(eid);
        if((*node_Gdim)[e.getIDs<Node>()[0]] == 3 && (*node_Gdim)[e.getIDs<Node>()[0]] == 3){
            (*edge_Gdim)[eid] = 3;
            cpt_out++;
        }
        else{
            //std::cout<<"in edge"<<std::endl;
            cpt_in++;
            (*edge_Gdim)[eid] = 0;
        }

        if(e.get<Face>().size() == 3){
            std::cout<<"Edge 3 faces"<<std::endl;
            exit(1);
        }
    }

    //std::cout<<"out = "<<cpt_out<<", in = "<<cpt_in<<std::endl;

    //Loop bord interne
    for(auto eid : m_smesh.edges()){
        Edge first_e = m_smesh.get<Edge>(eid);
        if((*edge_Gdim)[eid] == 0 && first_e.get<Face>().size() == 1){
            std::vector<Edge> in_loop;
            std::vector<Node> loop_nodes;
            std::vector<math::Point> points;

            //std::cout<<"On a trouvé un trou a "<<eid<<std::endl;
            //std::cout<<"A coté de la face "<<first_e.getIDs<Face>()[0]<<std::endl;
            bool closed = false;
            bool next = true;
            Edge current_e = first_e;
            int prev_n = current_e.get<Node>()[0].id();
            int cpt_it = 0;

            int cpt_manifold = 0;
            int n_manifold = -1;
            int i_manifold = -1;

            while(next){
                //std::cout<<"nb it "<<cpt_it++<<std::endl;
                cpt_it++;
                next = false;
                closed = false;
                in_loop.push_back(current_e);
                Node n = current_e.getIDs<Node>()[0] != prev_n ? current_e.get<Node>()[0]: current_e.get<Node>()[1];
                loop_nodes.push_back(n);
                points.push_back(n.point());
                for(auto n_e : n.get<Edge>()){
                    if((*edge_Gdim)[n_e.id()] == 0 && n_e.get<Face>().size() == 1 && n_e.id() != current_e.id()){
                        //if(!closed){
                        closed = true;
                        next = true;
                        prev_n = n.id();
                        current_e = n_e;
                        break;
                        //}
                        cpt_manifold++;
                    }
                }
                /*if(cpt_manifold > 1){
                    n_manifold = n.id();
                    i_manifold = cpt_it-1;
                }*/
                if(current_e.id() == eid){
                    break;
                }
                if(cpt_it > 10){
                    next = false;
                }
            }
            //std::cout<<"Loop terminée de taille "<<in_loop.size()<<std::endl;
            if(in_loop.size() > 2){

                math::Point center = math::Point::massCenter(points);
                Node center_n = m_smesh.newNode(center);
                for(auto e : in_loop){
                    Face tri = m_smesh.newTriangle(e.get<Node>()[0],e.get<Node>()[1],center_n);
                    e.add<Face>(tri);
                    //exit(0);
                }
            }else if(in_loop.size() == 2){
                Face tri = m_smesh.newTriangle(first_e.get<Node>()[0],first_e.get<Node>()[1],loop_nodes[1]);
                first_e.add<Face>(tri);
                in_loop[1].add<Face>(tri);
            }
            //break;
        }
    }

    //std::cout<<"Nb points on surf "<<m_smesh.getNbNodes()<<std::endl;
    //std::cout<<"Nb faces on surf "<<m_smesh.getNbFaces()<<std::endl;

    std::cout<<"Plan surface finished"<<std::endl;
    return &m_smesh;
}
/*----------------------------------------------------------------------------*/
