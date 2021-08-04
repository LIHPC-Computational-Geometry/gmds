//
// Created by calderans on 10/03/20.
//

#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/cad/FACSurface.h>


#include <gmds/graph/TopologyGraph.h>
#include <fstream>
#include <gmds/math/Orientation.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/

using namespace gmds;
using namespace graph;
/*----------------------------------------------------------------------------*/
TopologyGraph::TopologyGraph(cad::FACManager* Amanager, cad::GeomMeshLinker* Alinker, Mesh* Amesh):m_manager(Amanager),m_linker(Alinker),m_mesh(Amesh){
    BND_COLOR = m_mesh->getVariable<int,GMDS_FACE>("BND_SURFACE_COLOR");
    try {
        m_surface_color = m_mesh->newVariable<int, GMDS_FACE>("surface_color");
    }catch (GMDSException e) {
        m_surface_color = m_mesh->getVariable<int, GMDS_FACE>("surface_color");
    }
    try {
        m_isomorph = m_mesh->newVariable<int, GMDS_FACE>("isomorph");
    }catch (GMDSException e) {
        m_isomorph = m_mesh->getVariable<int, GMDS_FACE>("isomorph");
    }
    m_isomorph_tet = m_mesh->newVariable<int, GMDS_REGION>("isomorph_tet");
    for(auto f : m_mesh->faces()){
        (*m_surface_color)[f] = 0;
    }
    Variable<math::Vector3d> *var_X_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_X_POS");
    Variable<math::Vector3d> *var_Y_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_Y_POS");
    Variable<math::Vector3d> *var_Z_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_Z_POS");

    m_vertex_chart = m_mesh->getOrCreateVariable<math::Chart,GMDS_NODE>( "chart_field");

    m_face_weight = m_mesh->newVariable<double, GMDS_FACE>("face_weight");

    for(auto n : m_mesh->nodes()){
        math::Vector3d x = (*var_X_field)[n];
        math::Vector3d y = (*var_Y_field)[n];
        math::Vector3d z = (*var_Z_field)[n];
        math::Chart c(x,y,z);
        (*m_vertex_chart)[n]= c;
    }

    for(auto r:m_mesh->regions()){
        Region tet = m_mesh->get<Region>(r);
        for(auto f : tet.getIDs<Face>()){
            if(BND_COLOR->value(f) != 0){
                m_isomorph_tet->set(r,1);
            }
        }
    }
    m_cut_t = m_mesh->newVariable<int,GMDS_REGION>("CUT_tet");

    mark_face_surf = m_mesh->newMark<Face>();
    mark_tet_offset = m_mesh->newMark<Region>();

    m_cut_f = m_mesh->newVariable<int,GMDS_FACE>("CUT_face");

    m_minCut = new MinCut(Amesh);

    cut_id = 1;


}
/*----------------------------------------------------------------------------*/
int TopologyGraph::getGeomSurfID(TCellID AID){
    //Pour le moment on suppose qu'on pick un tet mais le mieux serait de pick un triangle
    int surf_ID = -1;
    Region r = m_mesh->get<Region>(AID);
    for(auto f : r.get<Face>()){
        if(f.get<Region>().size() == 1){
            //on ne sait pas si les faces sont classifiées donc on cherche la surface geom via les sommets
            /*for(auto n : f.get<Node>()){
                int cpt_surf_node = 0;
                if(m_linker->getGeomDim(n) == 3){
                    cpt_surf_node++;
                    if(cpt_surf_node == 3){
                        surf_ID = m_linker->getGeomId(n);
                        return surf_ID;
                    }
                }
            }*/
            if((*BND_COLOR)[f.id()] != 0) surf_ID = (*BND_COLOR)[f.id()];
        }
    }
    return surf_ID;
}
/*----------------------------------------------------------------------------*/
int TopologyGraph::generateCut(){
    if(!m_selection.empty()) {
        findSelectionBoundary();
        std::vector<TCellID> vide;
        boundaryMarking();

        IGMeshIOService ioService_m(m_mesh);
        VTKWriter vtkWriter_m(&ioService_m);

        computeWeight(vide);
        getCut();
        cutOptim();
        for(auto surf : m_selection){
            for(auto f : m_mesh->faces()){
                if(BND_COLOR->value(f) == surf){
                    m_cut_f->set(f,cut_id);
                }
            }
        }

        surfaceGraph();

        for(auto r : m_mesh->regions()){
            m_isomorph_tet->set(r,1);
        }

        vtkWriter_m.setCellOptions(gmds::N|gmds::F);
        vtkWriter_m.setDataOptions(gmds::N|gmds::F);
        vtkWriter_m.write("/ccc/temp/cont001/ocre/calderans/Results_debug/result_cut.vtk");


        cut_id++;
        m_selection.clear();
        m_boundaryEdges.clear();

        return 0;
    }
    return -1;
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::computeWeight(std::vector<TCellID> ATet){

    if(ATet.empty()) {
        for (auto f : m_mesh->faces()) {
            Face face = m_mesh->get<Face>(f);
            std::vector<math::Vector3d> closest_chart_vectors;
            std::vector<math::Point> points;
            std::vector<Node> nodes = face.get<Node>();

            //For each chart, get the most aligned vector with AV
            closest_chart_vectors.resize(3);
            points.resize(3);
            for (int i = 0; i < nodes.size(); i++) {
                Node n = nodes[i];
                math::Chart c = (*m_vertex_chart)[n.id()];

                math::Vector3d closest_vec;
                double max = 0;
                for (int i = 0; i < 3; i++) {
                    math::Vector3d vec = c[i];
                    double product = vec.dot(face.normal());
                    product = fabs(product);
                    if (product > max) {
                        max = product;
                        closest_vec = vec;
                    }
                }
                if (closest_vec.dot(face.normal()) < 0)
                    closest_vec = -1 * closest_vec;


                closest_chart_vectors[i] = closest_vec;
                points[i] = n.getPoint();
            }

            math::Point point = face.center();

            TCoord coords[3];
            math::Point::computeBarycentric(points[0], points[1], points[2], point,
                                            coords[0], coords[1], coords[2]);

            math::Vector3d norm = coords[0] * closest_chart_vectors[0];

            for (int i = 1; i < 3; i++) {
                norm += coords[i] * closest_chart_vectors[i];
            }
            norm.normalize();
            math::Vector3d norm2(1, 0, 0);

            double weight = 1 - fabs(norm.dot(face.normal()));
            if (weight > 0.90) { weight = 50; }
            else if (weight > 0.50) { weight = 12.5; }
            else if (weight > 0.25) { weight = 1; }

            (*m_face_weight)[f] = 1 + weight;

        }
    }else {
        for(auto t : ATet) {
            Region r = m_mesh->get<Region>(t);
            for (auto f : r.getIDs<Face>()) {
                Face face = m_mesh->get<Face>(f);
                std::vector<math::Vector3d> closest_chart_vectors;
                std::vector<math::Point> points;
                std::vector<Node> nodes = face.get<Node>();

                //For each chart, get the most aligned vector with AV
                closest_chart_vectors.resize(3);
                points.resize(3);
                for (int i = 0; i < nodes.size(); i++) {
                    Node n = nodes[i];
                    math::Chart c = (*m_vertex_chart)[n.id()];

                    math::Vector3d closest_vec;
                    double max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = c[i];
                        double product = vec.dot(face.normal());
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec = vec;
                        }
                    }
                    if (closest_vec.dot(face.normal()) < 0)
                        closest_vec = -1 * closest_vec;


                    closest_chart_vectors[i] = closest_vec;
                    points[i] = n.getPoint();
                }

                math::Point point = face.center();

                TCoord coords[3];
                math::Point::computeBarycentric(points[0], points[1], points[2], point,
                                                coords[0], coords[1], coords[2]);

                math::Vector3d norm = coords[0] * closest_chart_vectors[0];

                for (int i = 1; i < 3; i++) {
                    norm += coords[i] * closest_chart_vectors[i];
                }
                norm.normalize();
                math::Vector3d norm2(0, 0, 1);

                double weight = 1 - fabs(norm.dot(face.normal()));
                if (weight > 0.90) { weight = 100; }
                else if (weight > 0.50) { weight = 25; }
                else if (weight > 0.25) { weight = 2; }
                double weight2 =fabs(norm2.dot(face.normal()));
                if (weight2 > 0.9) weight2 = 1000;

                (*m_face_weight)[f] = 2 + weight + weight2;

            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::cutOptim(){

    int tet_treated = m_mesh->newMark<Region>();

    std::vector<TCellID> offset;
    std::vector<TCellID> faces;


    //profondeur 0
    for(auto f : m_mesh->faces()){
        Face face = m_mesh->get<Face>(f);
        if(face.get<Region>().size() == 2){
            if(m_cut_t->value(face.getIDs<Region>()[0]) != m_cut_t->value(face.getIDs<Region>()[1])){
                for(auto n : face.get<Node>()){
                    for(auto r : n.getIDs<Region>()){
                        if(!m_mesh->isMarked<Region>(r,tet_treated)){
                            offset.push_back(r);
                            m_mesh->mark<Region>(r,tet_treated);
                        }
                    }
                }
            }
        }
    }
    //profondeur 1
    std::vector<TCellID> new_offset;
    for(auto r : offset){
        Region region = m_mesh->get<Region>(r);
        for(auto n : region.get<Node>()){
            for(auto r : n.getIDs<Region>()){
                if(!m_mesh->isMarked<Region>(r,tet_treated)){
                    new_offset.push_back(r);
                    m_mesh->mark<Region>(r,tet_treated);
                }
            }
        }
    }
    for(auto t : new_offset){
        offset.push_back(t);
    }
    std::vector<TCellID> last_set;
    //profondeur 2
    for(auto r : new_offset){
        Region region = m_mesh->get<Region>(r);
        for(auto n : region.get<Node>()){
            for(auto r : n.getIDs<Region>()){
                if(!m_mesh->isMarked<Region>(r,tet_treated)){
                    m_mesh->mark<Region>(r,tet_treated);
                    last_set.push_back(r);
                }
            }
        }
    }

    //On prepare le cut
    for(auto t : last_set){
        offset.push_back(t);
        if(m_isomorph_tet->value(t) == 0) m_isomorph_tet->set(t,m_cut_t->value(t));
    }

    //On fait un truc moche pour le moment parce que sinon ça va être chiant
    int cpt_i = 0;
    std::map<TCellID,int> map;
    for(auto t : offset){
        m_mesh->mark<Region>(mark_tet_offset,t);
        map.emplace(t,cpt_i);
        cpt_i++;
    }

    computeWeight(offset);

    //Ici c'est pour recuperer les faces internes à l'offset
    std::set<TCellID> set_faces; //on utilise un set pour éviter les doublons
    for(auto r : offset){
        Region region = m_mesh->get<Region>(r);
        for(auto f : region.get<Face>()){
            if(f.getIDs<Region>().size() == 2){
                if(m_mesh->isMarked<Region>(f.getIDs<Region>()[0],tet_treated) && m_mesh->isMarked<Region>(f.getIDs<Region>()[1],tet_treated)){
                    set_faces.insert(f.id());
                }
            }
        }
    }
    //on transforme le set en vector
    for(auto f : set_faces){
        faces.push_back(f);
    }


    m_minCut->graph_cut(map,faces,m_isomorph_tet,m_cut_t,m_face_weight);

    m_mesh->unmarkAll<Region>(tet_treated);
    m_mesh->freeMark<Region>(tet_treated);
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> TopologyGraph::getGeomSurf(int AsurfID){
    cad::FACSurface* surf = m_manager->getFACSurface(AsurfID);
    std::set<TCellID> surf_tets;
    for(auto f: m_mesh->faces()) {
        if ((*BND_COLOR)[f] == surf->id()) {
            Face face = m_mesh->get<Face>(f);
            surf_tets.insert(face.getIDs<Region>()[0]);
        }
    }
    std::vector<TCellID> surf_vector;
    for(auto id: surf_tets){
        surf_vector.push_back(id);
    }
    return surf_vector;
}
/*----------------------------------------------------------------------------*/
int TopologyGraph::pickSurf(TCellID AID,std::vector<TCellID> &AIDsOutPut){
    int surf_id = getGeomSurfID(AID);
    if(surf_id == -1)
        return  -1;
    std::cout<<"surf id "<<surf_id<<std::endl;
    AIDsOutPut = getGeomSurf(surf_id);
    for(auto f : m_mesh->faces()){
        if((*BND_COLOR)[f] == surf_id){
            (*m_isomorph)[f] = 1;
        }
    }
    for(auto t:AIDsOutPut){
        (*m_isomorph_tet)[t] = 2;
    }
    m_selection.push_back(surf_id);

    return surf_id;
}
/*----------------------------------------------------------------------------*/
std::vector<int> TopologyGraph::findSelectionBoundary(){
    std::vector<int> selectionBoundary;
    for(auto s : m_selection){
        cad::FACSurface* surf = m_manager->getFACSurface(s);
        for(auto c : surf->curves()){
            int cpt = 0;
            for(auto s_c : c->surfaces()){
                int surf_id = s_c->id();
                if(std::find(m_selection.begin(), m_selection.end(), surf_id) != m_selection.end()){
                    cpt++;
                }
            }
            if(cpt == 1){
                selectionBoundary.push_back(c->id());
                std::vector<TCellID> mesh_edges;
                for(auto e : m_mesh->edges()){
                    Edge edge = m_mesh->get<Edge>(e);
                    if((m_linker->getGeomId(edge.get<Node>()[0]) == c->id() && m_linker->getGeomDim(edge.get<Node>()[0]) == 2 && ((m_linker->getGeomDim(edge.get<Node>()[1]) == 2 && m_linker->getGeomId(edge.get<Node>()[1]) == c->id()) || m_linker->getGeomDim(edge.get<Node>()[1]) == 1)) ||
                       (m_linker->getGeomId(edge.get<Node>()[1]) == c->id() && m_linker->getGeomDim(edge.get<Node>()[1]) == 2 && ((m_linker->getGeomDim(edge.get<Node>()[0]) == 2 && m_linker->getGeomId(edge.get<Node>()[0]) == c->id()) || m_linker->getGeomDim(edge.get<Node>()[0]) == 1)) ){
                        mesh_edges.push_back(e);
                    }
                }
                m_boundaryEdges.emplace(c->id(),mesh_edges);
            }
        }
    }
    std::cout<<"border size "<<selectionBoundary.size()<<std::endl;
    return selectionBoundary;
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::getCut(){

    m_minCut->graph_cut(m_isomorph_tet,m_cut_t,m_face_weight);
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::surfaceReconstruction(std::vector<TCellID> APath) {

    //Ici boucle recréation de la surface
    for(auto origin : m_mesh->faces()){
        //On cherche une face qui appartient à la coupe
        if(m_isomorph->value(origin) == 1) {
            std::set<int> adj_geom_surfaces;                    //Les faces géométriques adj à la nouvelle face
            std::vector<TCellID> new_surface_faces;             //L'ensemble des faces triangulaire d'une des nouvelles surfaces
            std::set<TCellID> next_faces_set;
            std::set<int> adj_geom_edges;

            next_faces_set.insert(origin);

            while (!next_faces_set.empty()) {

                std::vector<TCellID> next_faces_vec;
                next_faces_vec.reserve(next_faces_set.size());
                for (auto f : next_faces_set) {
                    next_faces_vec.push_back(f);
                }
                next_faces_set.clear();

                for (auto f_i : next_faces_vec) {

                    Face face = m_mesh->get<Face>(f_i);
                    new_surface_faces.push_back(f_i);
                    m_isomorph->set(f_i, 0);
                    m_mesh->mark<Face>(f_i, mark_face_surf);

                    for (auto e : face.get<Edge>()) {
                        if (std::find(APath.begin(), APath.end(), e.id()) ==
                            APath.end()) { //Si l'arête n'appartient pas à la nouvelle arête geom on continue l'algo
                            TCellID n0 = e.getIDs<Node>()[0], n1 = e.getIDs<Node>()[1];

                            //On teste les sommets de l'arete du maillage pour savoir si l'arêtes est sur une courbe geom
                            //Dans le future on va classifier les aretes du maillage directement pour gagner du temps de traitement
                            if (m_linker->getGeomDim<Node>(n0) == 2) {
                                int curve_id = m_linker->getGeomId<Node>(n0);
                                adj_geom_edges.insert(curve_id);
                            } else if (m_linker->getGeomDim<Node>(n1) == 2) {
                                int curve_id = m_linker->getGeomId<Node>(n1);
                                adj_geom_edges.insert(curve_id);
                            }


                            for (auto f : e.getIDs<Face>()) {
                                if (m_isomorph->value(f) != 0) {
                                    next_faces_set.insert(f);
                                }
                            }

                        }
                    }
                }
            }


            for (auto e : adj_geom_edges) {
                std::vector<TCellID> edges = m_boundaryEdges[e];
                int curvature = geomEdgeReclassification(
                        edges); //Normalement on obtient la classification de la courbe sur la nouvelle geométrie

                if (curvature ==
                    2) { //Si c'est flat (2) alors on cherche la surface adj qui ne fait pas parti de la selection
                    cad::GeomCurve *c = m_manager->getCurve(e);
                    for (auto s : c->surfaces()) {
                        int surf_id = s->id();
                        if (std::find(m_selection.begin(), m_selection.end(), surf_id) ==
                            m_selection.end()) {
                            adj_geom_surfaces.insert(surf_id);
                        }
                    }
                }
            }
            std::vector<int> adj_geom_surfaces_vec;
            adj_geom_surfaces_vec.reserve(adj_geom_surfaces.size());
            for (auto s : adj_geom_surfaces) {
                adj_geom_surfaces_vec.push_back(s);
            }
            //Si on a une seule surface flat adj alors on ajoute à cette surface les nouvelles faces
            if (adj_geom_surfaces_vec.size() == 1) {
                for (auto f : new_surface_faces) {
                    BND_COLOR->set(f, adj_geom_surfaces_vec[0]);
                    Face face = m_mesh->get<Face>(f);
                    for (auto n : face.getIDs<Node>()) {
                        int dim = m_linker->getGeomDim<Node>(n);
                        if (dim != 1 && dim != 2) {
                            m_linker->linkNodeToSurface(n, adj_geom_surfaces_vec[0]);
                        }
                    }
                }
            } else if (adj_geom_surfaces_vec.size() > 1) {
                int min = 999999;
                for (auto s : adj_geom_surfaces_vec) {
                    if (s < min) {
                        min = s;
                    }
                }
                for (auto f : new_surface_faces) {
                    BND_COLOR->set(f, min);
                    Face face = m_mesh->get<Face>(f);
                    for (auto n : face.getIDs<Node>()) {
                        int dim = m_linker->getGeomDim<Node>(n);
                        if (dim != 1 && dim != 2) {
                            m_linker->linkNodeToSurface(n, min);
                        }
                    }
                }
                for (auto f : m_mesh->faces()) {
                    if (BND_COLOR->value(f) != min &&
                        std::find(adj_geom_surfaces_vec.begin(), adj_geom_surfaces_vec.end(),
                                  BND_COLOR->value(f)) != adj_geom_surfaces_vec.end()) {
                        BND_COLOR->set(f, min);
                        Face face = m_mesh->get<Face>(f);
                        for (auto n : face.getIDs<Node>()) {
                            int dim = m_linker->getGeomDim<Node>(n);
                            if (dim == 3) {
                                m_linker->linkNodeToSurface(n, min);
                            }
                        }
                    }
                }
            } else {
                std::cout << "ERROR NOT IMPLEMENTED YET" << std::endl;
                std::cout << "\tNo adjacent flat surface found" << std::endl;
            }

        }
    }
}
/*----------------------------------------------------------------------------*/
int TopologyGraph::geomEdgeReclassification(std::vector<TCellID> AEdges){
    int nb_flat =0;
    int nb_convex=0;
    int nb_concave=0;
    bool face_not_found = true;
    for(auto e_id:AEdges){
        Edge ei = m_mesh->get<Edge>(e_id);
        std::vector<Face> adj_faces = ei.get<Face>();
        std::pair<Face,Face> surf_faces;
        for(auto adj : adj_faces){
            if(BND_COLOR->value(adj.id()) != 0 && std::find(m_selection.begin(),m_selection.end(),BND_COLOR->value(adj.id())) == m_selection.end()){
                surf_faces.first = adj;
                face_not_found = false;
            }else if(m_mesh->isMarked<Face>(adj.id(),mark_face_surf)){
                surf_faces.second = adj;
                face_not_found = false;
            }
        }
        if(face_not_found) std::cout<<"ERREUR GEOM EDGE RECLASSIFICATION"<<std::endl;
        std::vector<Node> ns0 = surf_faces.first.get<Node>();
        math::Vector3d v0 = surf_faces.first.normal();
        math::Vector3d v1 = surf_faces.second.normal();
        math::Point c1  = surf_faces.second.center();
        bool convex=true;
        if(math::Orientation::orient3d(c1,ns0[0].getPoint(),
                                       ns0[1].getPoint(),
                                       ns0[2].getPoint())==math::Orientation::NEGATIVE){
            //its non-convex
            convex=false;
        }
        double a = math::Constants::INVPIDIV180*v0.angle(v1);
        if(v0.dot(v1) < 0){
            a = fabs(a-180);
        }

        if(!convex)
            a+=180;
        if(a<45)
            nb_flat++;
        else if(60<a && a<120)
            nb_convex++;
        else if( a> 220)
            nb_concave++;
    }

    if(nb_convex>nb_concave && nb_convex>nb_flat)
        return 0;
    if(nb_concave>nb_convex && nb_concave>nb_flat)
        return 1;
    return 2;
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::computeEdgeFrameWeight(Edge e, math::Vector3d AtargetVec) {
    math::Vector3d edge_vec(e.get<Node>()[0].getPoint(),e.get<Node>()[1].getPoint());
    edge_vec.normalize();

    std::vector<math::Vector3d> closest_chart_vectors;
    std::vector<math::Point> points;
    std::vector<Node> nodes = e.get<Node>();

//For each chart, get the most aligned vector with AV
    closest_chart_vectors.resize(2);
    points.resize(2);
    for(int i = 0; i<nodes.size(); i++){
        Node n = nodes[i];
        math::Chart c = (*m_vertex_chart)[n.id()];

        math::Vector3d closest_vec;
        double max  = 0;
        for(int i = 0;i < 3; i++)
        {
            math::Vector3d vec = c[i];
            double product = vec.dot(AtargetVec);
            product = fabs(product);
            if(product > max)
            {
                max = product;
                closest_vec = vec;
            }
        }
        if( closest_vec.dot(AtargetVec) < 0 )
            closest_vec = -1*closest_vec;


        closest_chart_vectors[i] = closest_vec;
        points[i] = n.getPoint();
    }

    math::Vector3d xp(points[0],points[1]);
    math::Vector3d v(points[0],points[1]);

    double alpha = xp.norm()/v.norm();

    math::Vector3d norm = (1-alpha)*closest_chart_vectors[0];
    norm +=alpha*closest_chart_vectors[1];
    norm.normalize();
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::surfaceGraph() {

    std::set<TCellID> nodes_for_graph_set;
    std::vector<TCellID> nodes_for_graph;
    std::map<TCellID, double> edge_weight;
    std::vector<TCellID> target_nodes;
    std::set<TCellID> target_nodes_set;


    for(auto f : m_mesh->faces()){
        m_isomorph->set(f,0);
        Face face = m_mesh->get<Face>(f);



        //On récupère les faces à l'interface de la coupe
        if(face.get<Region>().size() == 2){
            if((m_cut_t->value(face.getIDs<Region>()[0]) == 1 && m_cut_t->value(face.getIDs<Region>()[1]) == 2) ||
               (m_cut_t->value(face.getIDs<Region>()[0]) == 2 && m_cut_t->value(face.getIDs<Region>()[1]) == 1)){
                m_isomorph->set(f,1);
                for(auto e : face.get<Edge>()){
                    math::Vector3d edge_vec(e.get<Node>()[0].getPoint(),e.get<Node>()[1].getPoint());
                    edge_vec.normalize();

                    std::vector<math::Vector3d> closest_chart_vectors;
                    std::vector<math::Point> points;
                    std::vector<Node> nodes = e.get<Node>();

                    //For each chart, get the most aligned vector with AV
                    closest_chart_vectors.resize(2);
                    points.resize(2);
                    for(int i = 0; i<nodes.size(); i++){
                        Node n = nodes[i];
                        math::Chart c = (*m_vertex_chart)[n.id()];

                        math::Vector3d closest_vec;
                        double max  = 0;
                        for(int i = 0;i < 3; i++)
                        {
                            math::Vector3d vec = c[i];
                            double product = vec.dot(face.normal());
                            product = fabs(product);
                            if(product > max)
                            {
                                max = product;
                                closest_vec = vec;
                            }
                        }
                        if( closest_vec.dot(face.normal()) < 0 )
                            closest_vec = -1*closest_vec;


                        closest_chart_vectors[i] = closest_vec;
                        points[i] = n.getPoint();
                    }

                    math::Vector3d xp(points[0],points[1]);
                    math::Vector3d v(points[0],points[1]);

                    double alpha = xp.norm()/v.norm();

                    math::Vector3d norm = (1-alpha)*closest_chart_vectors[0];
                    norm +=alpha*closest_chart_vectors[1];
                    norm.normalize();

                    if(m_linker->getGeomDim(e.get<Node>()[0]) != 2){
                        if(m_linker->getGeomDim(e.get<Node>()[1]) != 1 || m_linker->getGeomDim(e.get<Node>()[1]) != 2){
                            edge_weight.emplace(e.id(),1);
                            nodes_for_graph_set.insert(e.get<Node>()[0].id());
                            nodes_for_graph_set.insert(e.get<Node>()[1].id());
                            if(m_linker->getGeomDim(e.get<Node>()[0]) == 1){
                                target_nodes_set.insert(e.get<Node>()[0].id());
                            }
                            if(m_linker->getGeomDim(e.get<Node>()[1]) == 1){
                                target_nodes_set.insert(e.get<Node>()[1].id());
                            }
                        }
                    }else if(m_linker->getGeomDim(e.get<Node>()[1]) != 2){
                        if(m_linker->getGeomDim(e.get<Node>()[0]) != 1 || m_linker->getGeomDim(e.get<Node>()[0]) != 2){
                            edge_weight.emplace(e.id(),1);
                            nodes_for_graph_set.insert(e.get<Node>()[0].id());
                            nodes_for_graph_set.insert(e.get<Node>()[1].id());
                            if(m_linker->getGeomDim(e.get<Node>()[0]) == 1){
                                target_nodes_set.insert(e.get<Node>()[0].id());
                            }
                            if(m_linker->getGeomDim(e.get<Node>()[1]) == 1){
                                target_nodes_set.insert(e.get<Node>()[1].id());
                            }
                        }
                    }else{
                        continue;
                    }
                }
            }
        }
    }
    for(auto n : nodes_for_graph_set){
        nodes_for_graph.push_back(n);
    }


    std::vector<int> sources;
    std::vector<TCellID> path_edges;

    for(auto n : target_nodes_set){
        cad::GeomPoint* p = m_manager->getPoint(m_linker->getGeomId<Node>(n));
        int nb_curves_out = 0;
        int nb_curves_in = 0;
        for(auto c : p->curves()){
            bool curve_out = true;
            for(auto s : c->surfaces()){
                if(std::find(m_selection.begin(),m_selection.end(),s->id()) != m_selection.end()){
                    curve_out = false;
                }
            }
            if(curve_out){nb_curves_out++;}
            else{nb_curves_in++;}
        }
        if(nb_curves_out >= 1 && nb_curves_in < 3){
            sources.push_back(n);
        }
    }

    while(!sources.empty()) {

        TCellID source = sources.back();
        sources.pop_back();

        std::vector<int> geom_nodes_to_erase;
        std::vector<cad::GeomCurve *> curves;
        m_manager->getCurves(curves);
        for (auto c:curves) {

            cad::GeomPoint *p0 = c->points()[0];
            int id0 = p0->id();

            cad::GeomPoint *p1 = c->points()[1];
            int id1 = p1->id();

            if (id0 == source) {
                geom_nodes_to_erase.push_back(id1);
            } else if (id1 == source) {
                geom_nodes_to_erase.push_back(id0);
            }
        }

        for (auto n : target_nodes_set) {
            target_nodes.push_back(n);
        }

        target_nodes.erase(std::find(target_nodes.begin(), target_nodes.end(), source));
        for (auto gn : geom_nodes_to_erase) {
            //std::cout <<"Geom node to erase "<< gn+1 << std::endl;
            int to_erase = -1;
            for (auto n : target_nodes) {
                if (m_linker->getGeomId<Node>(n) == gn) {
                    //std::cout<<"Node "<<n<<std::endl;
                    to_erase = n;
                }
            }
            if (to_erase != -1) {
                target_nodes.erase(std::find(target_nodes.begin(), target_nodes.end(), to_erase));
            }
        }

        std::vector<TCellID> path = m_minCut->shortest_path(source, nodes_for_graph, edge_weight, target_nodes);

        sources.erase(std::find(sources.begin(),sources.end(),path[0]));

        //cad::GeomCurve* c = m_manager->newCurve();
        //int c_id = c->id();

        //Comme on a que les noeuds du maillage dans le path on cherche maintenant le chemin d'aretes correspondant
        for (auto e : m_mesh->edges()) {
            Edge edge = m_mesh->get<Edge>(e);
            for (int i = 0; i < path.size() - 1; i++) {
                if ((edge.getIDs<Node>()[0] == path[i] && edge.getIDs<Node>()[1] == path[i + 1]) ||
                    (edge.getIDs<Node>()[1] == path[i] && edge.getIDs<Node>()[0] == path[i + 1])) {
                    //On link les noeuds sur la nouvelle arete geom s'ils ne sont pas des sommets géométriques
                    if(m_linker->getGeomDim<Node>(path[i]) != 1){
                        //m_linker->linkNodeToCurve(path[i], c_id);
                    }
                    if(m_linker->getGeomDim<Node>(path[i+1]) != 1){
                        //m_linker->linkNodeToCurve(path[i+1], c_id);
                    }
                    path_edges.push_back(e);
                    break;
                }
            }
        }
        target_nodes.clear();
    }


    surfaceReconstruction(path_edges);

    m_mesh->unmarkAll<Face>(mark_face_surf);

}
/*----------------------------------------------------------------------------*/
int TopologyGraph::undoLastCut(){

    for(auto r:m_mesh->regions()){
        if(m_cut_t->value(r) == cut_id){
            m_cut_t->set(r,0);
        }
    }
    for(auto f:m_mesh->faces()){
        if(m_cut_f->value(f) == cut_id){
            m_cut_f->set(f,0);
        }
    }
    cut_id = cut_id-1;
}
/*----------------------------------------------------------------------------*/
int TopologyGraph::resetCut(){
    for(auto r:m_mesh->regions()){
        if(m_cut_t->value(r) != 0){
            m_cut_t->set(r,0);
        }
    }
    for(auto f:m_mesh->faces()){
        if(m_cut_f->value(f) != 0){
            m_cut_f->set(f,0);
        }
    }
    cut_id = 0;
}
/*----------------------------------------------------------------------------*/
int TopologyGraph::undoSelection(){

    int surf_id = m_selection.back();
    m_selection.pop_back();

    std::vector<TCellID> tet_surf = getGeomSurf(surf_id);
    for(auto f : m_mesh->faces()){
        if(BND_COLOR->value(f) == surf_id){
            m_isomorph->set(f,1);
        }
    }
    for(auto t:tet_surf){
        m_isomorph_tet->set(t,1);
    }
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::boundaryMarking(){

    //On stocke les surfaces déjà traitées pour éviter les opérations inutiles
    std::vector<int> id_surf_done;
    //tous les tets sont à m_isomorph = 1, sauf la selection
    for(auto b_e : m_boundaryEdges){
        cad::GeomCurve *c = m_manager->getCurve(b_e.first);
        if(c->getCurvatureInfo() == cad::GeomCurve::CONVEX){
            for(auto s : c->surfaces()){
                //Si la surface est convex et pas dans la selection on va marquer les tetras et les faces à 0
                if(std::find(m_selection.begin(),m_selection.end(),s->id()) == m_selection.end() &&
                   std::find(id_surf_done.begin(),id_surf_done.end(),s->id()) == id_surf_done.end()){
                    id_surf_done.push_back(s->id());
                    std::vector<TCellID> surf_tets = getGeomSurf(s->id());
                    for(auto f : m_mesh->faces()){
                        if((*BND_COLOR)[f] == s->id()){
                            (*m_isomorph)[f] = 0;
                        }
                    }
                    for(auto t:surf_tets){
                        (*m_isomorph_tet)[t] = 0;
                    }
                }
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::cutDirectionWeight(int ASurfID){

    int mark_treated = m_mesh->newMark<Face>();

    //Map des faces et de leur directions de coupe
    std::map<TCellID,math::Vector3d> directions;
    std::vector<TCellID> last_faces;
    for(auto f : m_mesh->faces()){
        //Faces au bord, donc on récupere le seul tetra incident
        if(BND_COLOR->value(f) == ASurfID){
            Face face = m_mesh->get<Face>(f);
            directions.emplace(f, face.normal());
            m_mesh->mark<Region>(f,mark_treated);
            last_faces.push_back(f);
        }
    }
    bool newEntries = true;
    do{
        std::vector<TCellID> new_faces;
        newEntries = false;
        for(auto f : last_faces){
            Face face = m_mesh->get<Face>(f);
            for(auto r : face.get<Region>()){
                if(m_mesh->isMarked<Region>(r.id(), mark_tet_offset)){
                    for(auto r_f : r.get<Face>()){
                        if(!m_mesh->isMarked<Face>(r_f.id(), mark_treated)){
                            newEntries = true;
                            m_mesh->mark<Face>(r_f.id(), mark_treated);
                            math::Vector3d last_norm = directions[f];

                            std::vector<math::Vector3d> closest_chart_vectors;
                            std::vector<math::Point> points;
                            std::vector<Node> nodes = r_f.get<Node>();

                            closest_chart_vectors.resize(3);
                            points.resize(3);
                            for (int i = 0; i < nodes.size(); i++) {
                                Node n = nodes[i];
                                math::Chart c = (*m_vertex_chart)[n.id()];

                                math::Vector3d closest_vec;
                                double max = 0;
                                for (int i = 0; i < 3; i++) {
                                    math::Vector3d vec = c[i];
                                    double product = vec.dot(last_norm);
                                    product = fabs(product);
                                    if (product > max) {
                                        max = product;
                                        closest_vec = vec;
                                    }
                                }
                                if (closest_vec.dot(last_norm) < 0)
                                    closest_vec = -1 * closest_vec;


                                closest_chart_vectors[i] = closest_vec;
                                points[i] = n.getPoint();
                            }

                            math::Point point = face.center();

                            TCoord coords[3];
                            math::Point::computeBarycentric(points[0], points[1], points[2], point,
                                                            coords[0], coords[1], coords[2]);

                            math::Vector3d norm = coords[0] * closest_chart_vectors[0];

                            for (int i = 1; i < 3; i++) {
                                norm += coords[i] * closest_chart_vectors[i];
                            }
                            norm.normalize();

                            directions.emplace(r_f.id(), norm);

                            new_faces.push_back(r_f.id());
                        }
                    }
                }
            }
        }

        last_faces = new_faces;

    }while(!newEntries);
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::noCut() {

    for(auto t : m_mesh->regions()){
        m_cut_t->set(t,1);
    }
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::cutFace(TCellID ANode1, TCellID ANode2){

    std::cout<<"ID node 1 = "<<m_linker->getGeomId<Node>(ANode1)<<std::endl;
    std::cout<<"Dim node 1 = "<<m_linker->getGeomDim<Node>(ANode1)<<std::endl;
    std::cout<<"ID node 2 = "<<m_linker->getGeomId<Node>(ANode2)<<std::endl;
    std::cout<<"Dim node 2 = "<<m_linker->getGeomDim<Node>(ANode1)<<std::endl;

    std::vector<cad::GeomSurface*> surfs1;
    std::vector<cad::GeomSurface*> surfs2;


    // On récupère les entités GEOM
    if(m_linker->getGeomDim<Node>(ANode1) == 1){
        cad::GeomPoint* p1 = m_manager->getPoint(m_linker->getGeomId<Node>(ANode1));
        surfs1 = p1->surfaces();
    } else if (m_linker->getGeomDim<Node>(ANode1) == 2){
        cad::GeomCurve* c1 = m_manager->getCurve(m_linker->getGeomId<Node>(ANode1));
        surfs1 = c1->surfaces();
    }

    if(m_linker->getGeomDim<Node>(ANode2) == 1){
        cad::GeomPoint* p2 = m_manager->getPoint(m_linker->getGeomId<Node>(ANode2));
        surfs2 = p2->surfaces();
    } else if (m_linker->getGeomDim<Node>(ANode2) == 2){
        cad::GeomCurve* c2 = m_manager->getCurve(m_linker->getGeomId<Node>(ANode2));
        surfs2 = c2->surfaces();
    }

    // On récupère la surface GEOM commune
    cad::GeomSurface* s;

    for(auto s1 : surfs1){
        for(auto s2 : surfs2){
            if(s1 == s2){
                s = s1;
            }
        }
    }
    int surf_id = s->id();
    std::cout<<"Common face "<<surf_id<<std::endl;

    // L'espace de travail est l'ensemble des arêtes du maillage tet classifiées sur la surface
    std::set<TCellID> edge_set;
    std::set<TCellID> node_set;
    for(auto f : m_mesh->faces()){
        if(m_linker->getGeomId<Face>(f) == surf_id){
            Face face = m_mesh->get<Face>(f);
            for(auto e : face.getIDs<Edge>()){
                edge_set.insert(e);
            }
            for(auto n : face.getIDs<Node>()){
                node_set.insert(n);
            }
        }
    }
    std::map<TCellID, double> edges_weight;
    for(auto e : edge_set){
        edges_weight.emplace(e,1);
    }
    std::vector<TCellID> nodes;
    for(auto n : node_set){
        nodes.push_back(n);
    }
    std::vector<TCellID> targetNode;
    targetNode.push_back(ANode2);

    // On set les variables de poids


    // On lance le shortest path entre les deux
    std::vector<TCellID> path = m_minCut->shortest_path(ANode1,nodes,edges_weight,targetNode);
    std::vector<TCellID> path_edges;

    //Comme on a que les noeuds du maillage dans le path on cherche maintenant le chemin d'aretes correspondant
    for (auto e : m_mesh->edges()) {
        Edge edge = m_mesh->get<Edge>(e);
        for (int i = 0; i < path.size() - 1; i++) {
            if ((edge.getIDs<Node>()[0] == path[i] && edge.getIDs<Node>()[1] == path[i + 1]) ||
                (edge.getIDs<Node>()[1] == path[i] && edge.getIDs<Node>()[0] == path[i + 1])) {
                //On link les noeuds sur la nouvelle arete geom s'ils ne sont pas des sommets géométriques
                if(m_linker->getGeomDim<Node>(path[i]) != 1){
                    //m_linker->linkNodeToCurve(path[i], c_id);
                }
                if(m_linker->getGeomDim<Node>(path[i+1]) != 1){
                    //m_linker->linkNodeToCurve(path[i+1], c_id);
                }
                path_edges.push_back(e);
                break;
            }
        }
    }

    // On crée une nouvelle courbe GEOM
    cad::GeomCurve* c = m_manager->newCurve();
    int c_id = c->id();

    for(int i = 1; i<path.size(); i++){
        //On va de 1 à path.size-1 pour pas prendre les sommets aux extrémités
        Node node = m_mesh->get<Node>(path[i]);
        m_linker->linkToCurve(node,c_id);
    }



    // On scinde la surface GEOM en deux
    std::set<TCellID> faces_new_surf;
    for(auto f : m_mesh->faces()){
        if(m_linker->getGeomId<Face>(f) == surf_id && m_linker->getGeomDim<Face>(f) == 3){
            m_mesh->mark<Face>(f,mark_face_surf);
            faces_new_surf.insert(f);
            Face origin = m_mesh->get<Face>(f);
            std::set<TCellID> next_surfs;
            for(auto e : origin.get<Edge>()){
                if(std::find(path_edges.begin(),path_edges.end(),e.id()) == path_edges.end()) {
                    for (auto f_e : e.getIDs<Face>()) {
                        if (m_linker->getGeomId<Face>(f_e) == surf_id && m_linker->getGeomDim<Face>(f_e) == 3 &&
                            f_e != f) {

                            next_surfs.insert(f_e);
                        }
                    }
                }
            }
            std::set<TCellID> tmp_set;
            while(!next_surfs.empty()){
                for(auto w_f : next_surfs){

                    m_mesh->mark<Face>(w_f,mark_face_surf);
                    faces_new_surf.insert(w_f);
                    Face current = m_mesh->get<Face>(w_f);
                    for(auto e : current.get<Edge>()){
                        if(std::find(path_edges.begin(),path_edges.end(),e.id()) == path_edges.end()) {
                            for (auto f_e : e.getIDs<Face>()) {
                                if (m_linker->getGeomId<Face>(f_e) == surf_id && m_linker->getGeomDim<Face>(f_e) == 3 &&
                                    f_e != w_f && !m_mesh->isMarked<Face>(f_e, mark_face_surf)) {
                                    tmp_set.insert(f_e);
                                }
                            }
                        }
                    }
                }
                next_surfs = tmp_set;
                tmp_set.clear();
            }
            break;
        }
    }

    cad::GeomSurface* new_surf = m_manager->newSurface();

    int id_new_surf = new_surf->id();
    std::cout<<"id surf = "<<id_new_surf<<" nb surfs "<<m_manager->getNbSurfaces()<<std::endl;
    for(auto f : faces_new_surf){

        Face face = m_mesh->get<Face>(f);
        m_linker->linkToSurface(face,id_new_surf);
        BND_COLOR->set(f,id_new_surf);
    }

    // On scinde les courbes GEOM des deux sommets GEOM
    // On reclassifie tout
}














