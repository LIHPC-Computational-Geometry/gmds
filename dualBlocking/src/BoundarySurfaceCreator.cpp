//
// Created by calderans on 07/03/19.
//

#include "gmds/dualBlocking/BoundarySurfaceCreator.h"
#include <map>
#include <gmds/math/Plane.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>

using namespace gmds;

using namespace db;


/*----------------------------------------------------------------------------*/

BoundarySurfaceCreator::BoundarySurfaceCreator(Mesh *AMesh,
        //TCellID AID,
        //int ASheetID,
                                               gmds::Variable<std::vector<DualSheetCreator::intersectInfo>>* info_Var,
                                               gmds::Variable<gmds::math::Chart>* vertex_chart_Var,
                                               gmds::Variable<int>* propagation_round_Var,
                                               gmds::Variable<int>* singularities_Var,
                                               int wave_precedent_mark,
                                               int face_treated_mark,
                                               int wave_tet_mark,
                                               double AMin_length,
                                               cad::GeomMeshLinker &ALinker,
                                               cad::FACManager &AManager):DualSheetCreator(AMesh,/* AID, ASheetID,*/ info_Var, vertex_chart_Var, propagation_round_Var, singularities_Var, wave_precedent_mark,face_treated_mark,wave_tet_mark,AMin_length,ALinker,AManager){


    wave_id = m_smesh.newVariable<int, GMDS_FACE>("wave");
    edge_Gdim = m_smesh.newVariable<int, GMDS_EDGE>("edge_dim");
    node_Gdim = m_smesh.newVariable<int, GMDS_NODE>("node_dim");
}
/*----------------------------------------------------------------------------*/
BoundarySurfaceCreator::~BoundarySurfaceCreator(){

}
/*----------------------------------------------------------------------------*/

bool BoundarySurfaceCreator::execute() {

    for (auto r : m_mesh->regions()) {
        (*m_propagation_round)[r] = -1;
    }
    surface.clear();

    std::cout<<"Tet "<<tetID<<std::endl;
    std::vector<TCellID > nodes = m_mesh->get<Region>(tetID).getIDs<Node>();
    for(auto n : nodes){
        if(linker.getGeomDim<Node>(n) == 3){
            std::cout<<"Surface "<<linker.getGeomId<Node>(n)<< " found"<<std::endl;
            setSurfaceID(linker.getGeomId<Node>(n));
        }
    }

    for(auto r : m_mesh->regions()){

        Region tet = m_mesh->get<Region>(r);

        std::vector<Node> nodes = tet.get<Node>();
        for(auto const &n : nodes){
            if(linker.getGeomId(n) == surface_id && linker.getGeomDim(n) == 3){
                bool good = false;
                math::Vector3d v;

                std::vector<Face> faces = tet.get<Face>();
                for(auto const &f : faces){
                    if(linker.getGeomId(f) == surface_id && linker.getGeomDim(f) == 3){
                        v = f.normal();
                        good = true;
                    }
                }

                if(good) {
                    (*m_propagation_round)[r] = 0;
                    surface.emplace(r, v);

                    std::vector<Edge> edges = tet.get<Edge>();
                    for (auto const &e : edges) {
                        m_mesh->mark(e, wave_precedent);
                    }
                }
            }
        }
    }


    std::map<gmds::TCellID,gmds::math::Vector3d> surface_tmp;

    std::set<TCellID> wave_tet, next_tet;
    std::vector<TCellID> wave_vec;

    surface_tmp = surface;

    for(auto const &t : surface_tmp){
        std::vector<Face> faces = m_mesh->get<Region>(t.first).get<Face>();
        for(auto const &f : faces){
            if(linker.getGeomId(f) == surface_id && linker.getGeomDim(f) == 3){
                std::vector<Node> nodes = f.get<Node>();
                for(auto n: nodes){
                    std::vector<TCellID> tets = n.getIDs<Region>();
                    for(auto tet : tets){
                        if((*m_propagation_round)[tet] == -1){
                            surface.emplace(tet,t.second);
                            (*m_propagation_round)[tet] = 1;
                        }
                    }
                }
            }
        }
    }

    /*for(auto const &t : surface){

        std::vector<Node> nodes = m_mesh->get<Region>(t.first).get<Node>();
        int cpt_node_on_surf = 0;
        for(auto const &n : nodes){
            if(linker.getGeomId(n) == surface_id && linker.getGeomDim(n) == 3) cpt_node_on_surf++;

            if(cpt_node_on_surf == 2){
                (*m_propagation_round)[t.first] = 0;
                break;
            }

            if(linker.getGeomDim(n) == 2){
                (*m_propagation_round)[t.first] = 1;
            }
        }

        if((*m_propagation_round)[t.first] == 1){
            wave_tet.insert(t.first);
            wave_vec.push_back(t.first);
        }
    }

    std::vector<intersectInfo> infos;
    std::set<TCellID> nodes_moved;

    std::set<Edge> edges_set;

    //return true;

    while(!wave_tet.empty()) {
        for (auto t : wave_tet) {

            bool moved_node_init = false;
            TCellID moved_node_id_init = NullID;

            Region tet = m_mesh->get<Region>(t);
            math::Plane plane(tet.center(), surface[t]);

            std::vector<intersectInfo> infos_tmp = intersect(t, plane, moved_node_init, moved_node_id_init);
            if(moved_node_init){
                nodes_moved.insert(moved_node_id_init);
            }

            for (auto const &i : infos_tmp) {
                //std::cout << "i.e = " << i.e << std::endl;
                (*var_info)[i.e].push_back(i);
                m_mesh->mark(m_mesh->get<Edge>(i.e), wave_precedent);
            }
        }

        for (auto n : nodes_moved) {

            Node node_moved = m_mesh->get<Node>(n);

            std::vector<Region> regions_to_update = node_moved.get<Region>();

            for(auto const &r : regions_to_update){
                if((*m_propagation_round)[r.id()] == 2){
                    next_tet.insert(r.id());

                    std::vector<Edge> edges = r.get<Edge>();
                    for(auto const &e : edges){
                        (*var_info)[e.id()].clear();
                    }
                }else if((*m_propagation_round)[r.id()] == 0) continue;
            }

            std::vector<Edge> adj_edges = node_moved.get<Edge>();
            double d = 10000;
            for (auto const &a : adj_edges) {
                std::vector<Node> a_nodes = a.get<Node>();
                Node opp_node = (a_nodes[0].id() == n) ? a_nodes[1] : a_nodes[0];
                double ad = opp_node.getPoint().distance(node_moved.getPoint());
                if (ad < d)
                    d = ad;
            }

            d = d / (5 * sqrt(3));
            double dx = drand48() * d;
            double dy = drand48() * d;
            double dz = drand48() * d;
            //std::cout << " \t from: " << node_moved.getPoint() << std::endl;

            math::Point new_loc(node_moved.getPoint().X() + dx,
                                node_moved.getPoint().Y() + dy,
                                node_moved.getPoint().Z() + dz);
            //std::cout << " \t to new location: " << new_loc << std::endl;
            node_moved.setPoint(new_loc);
            //std::cout << "node moved = " << node_moved.id() << std::endl;
        }

        nodes_moved.clear();

        wave_tet = next_tet;
        next_tet.clear();
    }

    wave_tet.clear();
    std::set<intersectInfo> infos_set;

    for(auto const &t : wave_vec){

        (*m_propagation_round)[t] = 0;

        std::vector<Edge> edges = m_mesh->get<Region>(t).get<Edge>();
        for(auto const &e : edges){
            if(!(*var_info)[e.id()].empty()){
                infos_set.insert((*var_info)[e.id()][0]);
                std::vector<TCellID > tet_e = e.getIDs<Region>();
                for(auto const &t_e : tet_e){
                    if((*m_propagation_round)[t_e] == -1){
                        wave_tet.insert(t_e);
                        surface.emplace(t_e,(*var_info)[e.id()][0].v);
                        (*m_propagation_round)[t_e] = 0;
                    }
                }
            }
        }
    }*/

    //bool result = propagation_loop(wave_tet);

    return  true;
}

/*----------------------------------------------------------------------------*/

int BoundarySurfaceCreator::getSurfaceID(){
    return surface_id;
}
/*----------------------------------------------------------------------------*/

void BoundarySurfaceCreator::setSurfaceID(int AID){
    surface_id = AID;
}
/*----------------------------------------------------------------------------*/
Mesh* BoundarySurfaceCreator::buildSurfaceSheet() {

    m_smesh.clear();

    std::map<int,Node> node_mesh_to_smesh;
    std::set<Face> surface_faces;
int cpt_faces = 0;

    for(auto t : surface){
        Region r = m_mesh->get<Region>(t.first);
        std::vector<Face> faces = r.get<Face>();

        for(auto f : faces){
            cpt_faces++;
            if(f.getIDs<Region>().size() == 2) {
                TCellID ids0 = f.getIDs<Region>()[0];
                TCellID ids1 = f.getIDs<Region>()[1];

                if (((*m_propagation_round)[ids0] > -1 && (*m_propagation_round)[ids1] == -1) ||
                    ((*m_propagation_round)[ids1] > -1 && (*m_propagation_round)[ids0] == -1)) {
                    surface_faces.insert(f);
                }
            }
        }
    }

    std::set<TCellID> nodes_set;
    for(auto f : surface_faces){
        std::vector<TCellID> nodes = f.getIDs<Node>();
        for(auto n : nodes){
            nodes_set.insert(n);
        }
    }
    for(auto n : nodes_set){
        Node node = m_mesh->get<Node>(n);
        Node n_new = m_smesh.newNode(node.getPoint());
        node_mesh_to_smesh[n] = n_new;
            (*node_Gdim)[n_new.id()] = linker.getGeomDim<Node>(n);
    }
    for(auto f : surface_faces){
        std::vector<Node> nodes_face;
        std::vector<TCellID> nodes = f.getIDs<Node>();
        for(auto n : nodes){
            nodes_face.push_back(node_mesh_to_smesh[n]);
        }
        m_smesh.newFace(nodes_face);
    }

    return &m_smesh;
}
/*----------------------------------------------------------------------------*/


















