//
// Created by calderans on 10/07/19.
//

#include <Predicates_psm.h>
#include <gmds/math/Segment.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include "gmds/dualBlocking/DualSheetCreator.h"

using namespace gmds;
using namespace db;

DualSheetCreator::DualSheetCreator(Mesh *AMesh,
        // TCellID AID,
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
                                   cad::FACManager &AManager):m_mesh(AMesh),
                                   m_smesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                                   R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E)),
                                   linker(ALinker),manager(AManager) {

    //tetID = AID;
    //sheet_ID = ASheetID;

    propagation_round = 0;

    m_vertex_chart = vertex_chart_Var;

    var_info = info_Var;

    m_propagation_round = propagation_round_Var;

    min_length = AMin_length;

    m_singularities = singularities_Var;

    wave_precedent = wave_precedent_mark;
    face_treated = face_treated_mark;

    mark_wave_tet = wave_tet_mark;


}
/*----------------------------------------------------------------------------*/

DualSheetCreator::~DualSheetCreator(){

    //clear();
}
/*----------------------------------------------------------------------------*/

void DualSheetCreator::setTetID(gmds::TCellID AID){
    tetID = AID;
}
/*----------------------------------------------------------------------------*/

void DualSheetCreator::setSheetID(int AID){
    sheet_ID = AID;
}
/*----------------------------------------------------------------------------*/

int DualSheetCreator::getID(){
    return sheet_ID;
}
/*----------------------------------------------------------------------------*/

std::map<gmds::TCellID,gmds::math::Vector3d> DualSheetCreator::getSurface(){
    return surface;
}
/*----------------------------------------------------------------------------*/

bool DualSheetCreator::propagation_loop(std::set<TCellID> ATet_list) {


    IGMeshIOService m_ioService(m_mesh);
    Mesh imesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doci(&imesh);
    doci.buildFacesAndR2F();
    doci.buildEdgesAndX2E();
    doci.updateUpwardConnectivity();

    gmds::Variable<int>* m_blocks = m_mesh->getVariable<int, GMDS_REGION>("blocks");

    std::set<TCellID> wave_tet = ATet_list;

    std::set<TCellID> boundary_edges;

    for (auto t : wave_tet) {
        std::vector<Edge> tmp_edge = m_mesh->get<Region>(t).get<Edge>();
        for (auto const &e : tmp_edge) {

            std::vector<TCellID> tmp_tet = e.getIDs<Region>();
            for (auto tet : tmp_tet) {
                if ((*m_propagation_round)[tet] == -1) {
                    boundary_edges.insert(e.id());
                }
            }
        }
    }

    for (int i = 1; !wave_tet.empty(); i++) {

        std::vector<TCellID> wave_tet_vec;
        wave_tet_vec.reserve(wave_tet.size());
        for (auto t : wave_tet) {
            wave_tet_vec.push_back(t);
        }

        bool restart_wave_for_node_intersection;
        bool restart_wave_for_plane_issue;

        do {
            restart_wave_for_plane_issue = false;
            do {
                restart_wave_for_node_intersection = false;
                std::set<TCellID> nodes_to_move;

                for (auto i_tet = 0; i_tet < wave_tet_vec.size(); i_tet++) {


                    //To list the intersection points created in the tet for out put
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    std::vector<intersectInfo> for_output;
                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    TCellID t = wave_tet_vec[i_tet];


                    m_mesh->mark<Region>(t, mark_wave_tet);

                    std::vector<intersectInfo> from_cut;
                    std::vector<Edge> tmp_edge = m_mesh->get<Region>(t).get<Edge>();
                    for (auto const &e : tmp_edge) {
                        if (m_mesh->isMarked(e, wave_precedent)) {// we only make plan from previous wave edge
                            std::vector<intersectInfo> e_info = (*var_info)[e.id()];
                            if (!e_info.empty()) {
                                from_cut.insert(from_cut.end(), e_info.begin(), e_info.end());
                            }
                        }
                    }
                    //std::cout << "from cut size " << from_cut.size() << std::endl;

                    bool moved_node = false;
                    TCellID moved_node_id = NullID;


                    std::vector<intersectInfo> out_cut;
                    out_cut = face_cut(t, from_cut, moved_node, moved_node_id);

                    if (moved_node) {
                        //The point has not been moved, we only know that it should be to avoid issue
                        //But it must be done by taking care of all the tets around the edge we look at
                        nodes_to_move.insert(moved_node_id);
                        //std::cout<<"Node moved"<<std::endl;

                    } else {
                        //We update output edge infos
                        for (auto cut2 : out_cut) {
                            if (!m_mesh->isMarked<Edge>(cut2.e, wave_precedent)) {
                                //std::cout<<"Test 1 "<<2201<<":"<<linker.getGeomDim<Face>(2201)<<std::endl;
                                (*var_info)[cut2.e].push_back(cut2);
                                //std::cout<<"Test 2 "<<2201<<":"<<linker.getGeomDim<Face>(2201)<<std::endl;
                            }
                        }
                    }

                }
                wave_tet_vec.clear();
                //CORRECTION FOR NODE MOVEMENT
                for (auto n_id:nodes_to_move) {
                    //std::cout << "NODE TO MOVE: " << n_id << std::endl;

                    restart_wave_for_node_intersection = true;
                    //we restart only some tets
                    Node node_moved = m_mesh->get<Node>(n_id);
                    std::vector<Region> regions_to_update = node_moved.get<Region>();
                    int nb_tet_updated = 0;
                    for (auto const &r: regions_to_update) {
                        if (m_mesh->isMarked(r, mark_wave_tet)) {
                            m_mesh->unmark(r, mark_wave_tet);

                            //std::cout << "Region updated " << r.id() << std::endl;

                            nb_tet_updated++;
                            wave_tet_vec.push_back(r.id());
                        }
                    }
                    //std::cout << "Nb tet updated " << nb_tet_updated << std::endl;
                    std::vector<Face> faces_to_update = node_moved.get<Face>();
                    for (auto const &f : faces_to_update) {
                        if (m_mesh->isMarked(f, face_treated)) {
                            int nb_pre_edges = 0;
                            std::vector<Edge> f_edges = f.get<Edge>();
                            for (auto const &f_e: f_edges) {
                                if (m_mesh->isMarked(f_e, wave_precedent)) {
                                    nb_pre_edges++;
                                }
                            }
                            if (nb_pre_edges != 2) {

                                for (auto const &e_r : f_edges) {
                                    if (!m_mesh->isMarked(e_r, wave_precedent) && (e_r.getIDs<Node>()[0] == n_id
                                                                                   || e_r.getIDs<Node>()[1] == n_id)) {

                                        (*var_info)[e_r.id()].clear();
                                        m_mesh->unmark(f, face_treated);
                                    }
                                }
                            }
                        }
                    }
                    std::vector<Edge> adj_edges = node_moved.get<Edge>();
                    double d = 10000;
                    for (auto const &a: adj_edges) {
                        std::vector<Node> a_nodes = a.get<Node>();
                        Node opp_node = (a_nodes[0].id() == n_id) ? a_nodes[1] : a_nodes[0];
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
                }

            } while (restart_wave_for_node_intersection);

        }while (restart_wave_for_plane_issue);

        std::set<TCellID > to_erase;

        for(auto t : wave_tet){

            std::vector<math::Vector3d> vecs;

            std::vector<Edge> t_edges = m_mesh->get<Region>(t).get<Edge>();

            for (auto const &e : t_edges) {


                if (!(*var_info)[e.id()].empty()) {
                    //std::cout<<(*var_info)[e.getID()][0].wave<<std::endl;
                    if ((*var_info)[e.id()][0].wave >= propagation_round - 2) {
                        vecs.push_back((*var_info)[e.id()][0].v);
                    }
                }
            }

            if(!vecs.empty()) {

                math::Vector3d mean_vec = meanVector(vecs);

                for (auto const &edge : t_edges) {
                    std::vector<TCellID> neighbors = edge.getIDs<Region>();
                    for (auto n : neighbors) {
                        if ((*m_propagation_round)[n] < (*m_propagation_round)[t] - 2 &&
                            (*m_propagation_round)[n] != -1) {

                            //now we verify if the neighbour is "parallel" to the tet (used to auto intersecting sheets)
                            if(surface[n].dot(mean_vec) > 0.75) {
                                //std::cout<<"delete tetra "<<t<<" with round "<<(*m_propagation_round)[t]<<" adjacent to "<<n<<" with round "<<(*m_propagation_round)[n]<<std::endl;
                                to_erase.insert(t);
                            }
                        }
                    }
                }
            }
            else{
                //No new point in the tet
                to_erase.insert(t);
            }
        }


        for (auto t : wave_tet) {
            std::vector<Edge> tet_edge = m_mesh->get<Region>(t).get<Edge>();
            for (auto const &e : tet_edge) {
                if (!m_mesh->isMarked(e, wave_precedent)) {

                    //We do not consider input edges but only output and inner wave edges

                    std::vector<intersectInfo> e_info = (*var_info)[e.id()];
                    //if 1 point nothing to do
                    if (e_info.size() > 1) {

                        intersectInfo average_e_info;
                        averageInfo(e_info, average_e_info);
                        average_e_info.tet = e_info[0].tet;
                        math::Segment s(e.get<Node>()[0].getPoint(), e.get<Node>()[1].getPoint());

                        (*var_info)[e.id()].clear();
                        (*var_info)[e.id()].push_back(
                                average_e_info);//we set the intersection info to the average intersection info
                    }
                }
            }
        }


        for(auto t : wave_tet){

            math::Vector3d mean_vec;
            std::vector<math::Vector3d> vecs;

            std::vector<Edge> tmp_edge = m_mesh->get<Region>(t).get<Edge>();
            for (auto const &e : tmp_edge) {

                if (!(*var_info)[e.id()].empty()) {
                    if ((*var_info)[e.id()][0].wave >= propagation_round - 2) {
                        vecs.push_back((*var_info)[e.id()][0].v);
                    }
                }
            }

            if(vecs.empty()){std::cout<<"ERROR: not point in tet"<<std::endl;return false;}

            mean_vec = meanVector(vecs);

            surface.emplace(t,mean_vec);
        }

        //erase tet with issues so it doesnt impact the sheet later
        for(auto t : to_erase){
            wave_tet.erase(t);
        }

        boundary_edges.clear();


        for (auto t : wave_tet) {

            std::vector<Edge> tmp_edge = m_mesh->get<Region>(t).get<Edge>();
            for (auto const &e : tmp_edge) {

                std::vector<TCellID> tmp_tet = e.getIDs<Region>();
                for (auto tet : tmp_tet) {
                    if ((*m_propagation_round)[tet] == -1) {
                        boundary_edges.insert(e.id());
                    }else if(!(*var_info)[e.id()].empty()){


                        //for auto intersect
                        //----------------------------------------------
                        /*if((*var_info)[e.id()][0].v.dot(surface[tet]) < 0.70){
                            boundary_edges.insert(e.id());
                        }*/
                    }
                }
            }
        }

        wave_tet.clear();

        propagation_round++;

        /*for(auto element : surface){
            std::cout<<element.first<<" : propagation = "<<(*m_propagation_round)[element.first]<<", "<<element.second.X()<<", "<<element.second.Y()<<", "<<element.second.Z()<<std::endl;
        }*/


        for (auto e : boundary_edges) {

            Edge edge = m_mesh->get<Edge>(e);
            if ((*var_info)[e].size() == 1) {

                std::vector<TCellID> tmp_tet = edge.getIDs<Region>();
                for (auto t : tmp_tet) {

                    if ((*m_propagation_round)[t] == -1) {
                        if(!isAxisSingularity(t,(*var_info)[e][0].v)) {
                            wave_tet.insert(t);
                            (*m_propagation_round)[t] = propagation_round;
                            /* if(propagation_round == 8){
                              Node n1 = imesh.newNode(edge.get<Node>()[0].getPoint());
                              Node n2 = imesh.newNode(edge.get<Node>()[0].getPoint());
                              Node n3 = imesh.newNode(edge.get<Node>()[1].getPoint());

                              std::vector<Node> nodes;
                              nodes.push_back(n1);
                              nodes.push_back(n2);
                              nodes.push_back(n3);

                              imesh.newFace(nodes);
                              }*/
                        }
                        else{
                            return false;
                        }
                    }else if((*m_propagation_round)[t] != propagation_round){
                        //for auto intersect
                        //----------------------------------------------
                        /*
                        if((*var_info)[e][0].v.dot(surface[t]) < 0.70){

                            //std::cout<<"var info "<<std::endl;
                            //std::cout<<(*var_info)[e][0].v.X()<<", "<<(*var_info)[e][0].v.Y()<<", "<<(*var_info)[e][0].v.Z()<<std::endl;
                            std::cout<<"crossed in tet "<<t<<std::endl;
                            std::cout<<"propagation in tet "<<(*m_propagation_round)[t]<<std::endl;
                            std::cout<<surface[t].X()<<", "<<surface[t].Y()<<", "<<surface[t].Z()<<std::endl;


                            wave_tet.insert(t);
                            std::vector<TCellID> edges_id = m_mesh->get<Region>(t).getIDs<Edge>();
                            for(auto eid : edges_id){
                                if(!(*var_info)[eid].empty()){
                                    if((*var_info)[eid][0].v.dot(surface[t]) < 0.70){
                                        (*var_info)[eid].clear();
                                    }
                                }
                            }
                            std::vector<Face> faces= m_mesh->get<Region>(t).get<Face>();
                            for(auto f : faces){
                                m_mesh->unmark(f,face_treated);
                            }
                            (*m_propagation_round)[t] = propagation_round;
                        }*/
                    }
                }

                m_mesh->mark(edge, wave_precedent);

            } else if ((*var_info)[e].size() > 1) {
                Edge edge_fail = m_mesh->get<Edge>(e);
                std::cout<<"More than one point on edge "<<e<<" "<<edge_fail.getIDs<Node>()[0]<<","<<edge_fail.getIDs<Node>()[1]<<
                         " in tet "<<(*var_info)[e][0].tet <<std::endl;
                for(auto inf : (*var_info)[e]){
                    std::cout<<"info in tet "<<inf.tet<<", point ("<<inf.p.X()<<", "<<inf.p.Y()<<", "<<inf.p.Z()<<std::endl;
                }
                exit(12);
            }
        }



        m_mesh->unmarkAll<Region>(mark_wave_tet);

        /*gmds::VTKWriter m_vtkWriter(&m_ioService);
        m_vtkWriter.setCellOptions(gmds::N|gmds::R);
        m_vtkWriter.setDataOptions(gmds::N|gmds::R);
        m_vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/block_iterations"+std::to_string(i)+".vtk");*/
        //std::cout<<std::endl;
    }
    /*
    for(auto e : m_mesh->edges()){
      bool wave_8 = false;
      bool wave_7 = false;
      Edge edge = m_mesh->get<Edge>(e);
      std::vector<TCellID> colors = edge.getIDs<Region>();
      for(auto r : colors){
	if ((*m_propagation_round)[r] == 7) {
	  wave_7 = true;
	}
	if ((*m_propagation_round)[r] == 8) {
	  wave_8 = true;
	}
      }
      if(wave_7 && wave_8){
	Node n1 = imesh.newNode(edge.get<Node>()[0].getPoint());
	Node n2 = imesh.newNode(edge.get<Node>()[0].getPoint());
	Node n3 = imesh.newNode(edge.get<Node>()[1].getPoint());

	std::vector<Node> nodes;
	nodes.push_back(n1);
	nodes.push_back(n2);
	nodes.push_back(n3);

	imesh.newFace(nodes);
      }
      }*/

    //IGMeshIOService i_ioService(&imesh);
    //gmds::VTKWriter i_vtkWriter(&i_ioService);
    //i_vtkWriter.setCellOptions(gmds::N|gmds::F);
    //i_vtkWriter.setDataOptions(gmds::N|gmds::F);
    // i_vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/points.vtk");
    /*Mesh imesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doci(&imesh);
    doci.buildFacesAndR2F();
    doci.buildEdgesAndX2E();
    doci.updateUpwardConnectivity();

    for(auto tet : surface){
        std::vector<Edge> edges =m_mesh->get<Region>(tet.first).get<Edge>();
        for(auto e : edges){
            if(!(*var_info)[e.id()].empty()){
                Node n1 = imesh.newNode(e.get<Node>()[0].getPoint());
                Node n2 = imesh.newNode(e.get<Node>()[0].getPoint());
                Node n3 = imesh.newNode(e.get<Node>()[1].getPoint());

                std::vector<Node> nodes;
                nodes.push_back(n1);
                nodes.push_back(n2);
                nodes.push_back(n3);

                imesh.newFace(nodes);
            }
        }
    }

    IGMeshIOService ioService(&imesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/edges.vtk");*/



    return true;
}
/*----------------------------------------------------------------------------*/

const math::Vector3d DualSheetCreator::findNormal(math::Point& APoint, const math::Vector3d& AV, Region& ATet)
{
    //std::cout<<" ================= Finding normal ================= "<<std::endl;

    // We get the charts of the four nodes of ATet
    std::vector<math::Vector3d> closest_chart_vectors;
    std::vector<math::Point> points;
    std::vector<Node> nodes = ATet.get<Node>();

    //For each chart, get the most aligned vector with AV
    closest_chart_vectors.resize(4);
    points.resize(4);
    for(int i = 0; i<nodes.size(); i++){
        Node n = nodes[i];
        math::Chart c = (*m_vertex_chart)[n.id()];


        closest_chart_vectors[i] = closestComponentVector(AV, c);
        points[i] = n.getPoint();
    }
    TCoord coords[4];
    math::Point::computeBarycentric(points[0], points[1], points[2], points[3], APoint,
                                    coords[0], coords[1], coords[2],coords[3]);

    math::Vector3d norm = coords[0]*closest_chart_vectors[0];

    for(int i = 1; i<4; i++){
        norm +=coords[i]*closest_chart_vectors[i];
    }
    norm.normalize();


    return norm;
}
/*----------------------------------------------------------------------------*/

const math::Vector3d DualSheetCreator::findNormal(math::Point& APoint, const math::Vector3d& AV, Edge& AEdge)
{
    //std::cout<<" ================= Finding normal ================= "<<std::endl;

    // We get the charts of the two nodes of AEdge
    std::vector<math::Vector3d> closest_chart_vectors;
    std::vector<math::Point> points;
    std::vector<Node> nodes = AEdge.get<Node>();

    //For each chart, get the most aligned vector with AV
    closest_chart_vectors.resize(2);
    points.resize(2);
    for(int i = 0; i<nodes.size(); i++){
        Node n = nodes[i];
        math::Chart c = (*m_vertex_chart)[n.id()];
        closest_chart_vectors[i] = closestComponentVector(AV, c);
        points[i] = n.getPoint();
    }

    math::Vector3d xp(points[0],APoint);
    math::Vector3d v(points[0],points[1]);

    double alpha = xp.norm()/v.norm();

    math::Vector3d norm = (1-alpha)*closest_chart_vectors[0];
    norm +=alpha*closest_chart_vectors[1];
    norm.normalize();


    //std::cout<<"Normal = ("<<norm.X()<<","<<norm.Y()<<","<<norm.Z()<<")"<<std::endl;
    return norm;
}
/*----------------------------------------------------------------------------*/

const math::Vector3d DualSheetCreator::closestComponentVector(const math::Vector3d& AV,
                                                              math::Chart& AChart)
{
    //std::cout<<" ================= Finding closest component vector ================= "<<std::endl;

    math::Vector3d closest_vec;
    double max  = 0;
    for(int i = 0;i < 3; i++)
    {
        math::Vector3d vec = AChart[i];
        double product = vec.dot(AV);
        product = fabs(product);
        if(product > max)
        {
            max = product;
            closest_vec = vec;
        }
    }
    if( closest_vec.dot(AV) < 0 )
        closest_vec = -1*closest_vec;
    return closest_vec;
}

/*----------------------------------------------------------------------------*/

std::vector<DualSheetCreator::intersectInfo>
DualSheetCreator::intersect(TCellID ATetra, math::Plane APlane,
                            bool& AMoveNode,
                            gmds::TCellID& AMovedNodeID)
{
    //std::cout<<"================ intersecting ================"<<std::endl;
    Region tet = m_mesh->get<Region>(ATetra);

    std::vector<intersectInfo> infos;
    std::vector<Edge> edges = tet.get<Edge>();

    for (int i = 0; i < 6; i++) {

        if(!AMovedNodeID == 9999){
            std::cout<<"Continue"<<std::endl;
            if(!(*var_info)[edges[i].id()].empty()) continue;
        }
        intersectInfo iI;
        iI.wave = propagation_round;
        iI.tet = ATetra;
        math::Point PI;
        math::Vector3d v;
        iI.e = edges[i].id();
        int intersection = intersectEdge(APlane, PI, edges[i]);
        //std::cout<<"Intersection = "<<intersection<<std::endl;
        if ( intersection == 1) {

            iI.p = PI;
            iI.v = findNormal(PI, APlane.getNormal(), edges[i]);
            math::Plane plane(iI.p, iI.v);

            math::Point proj = plane.project(APlane.getPoint());
            math::Vector3d vec(proj, iI.p);
            iI.d = vec;

            infos.push_back(iI);

        } else if(intersection == -1){
            AMoveNode = true;
            AMovedNodeID = edges[i].get<Node>()[0].id();
            return infos;
        } else if(intersection == -2){
            std::cout<<"tetra "<<ATetra<<std::endl;
            AMoveNode = true;
            AMovedNodeID = edges[i].get<Node>()[1].id();
            return infos;
        }
        else{
            continue;
        }
    }

    return infos;
}
/*----------------------------------------------------------------------------*/

int DualSheetCreator::intersectEdge(math::Plane& APl, math::Point& PI,
                                    Edge& AE)
{
    //std::cout<<" ================= Intersecting Plane/Edge ================= "<<std::endl;
    //By default, we do not expect to move a node of the tet

    std::vector<Node> nodes = AE.get<Node>();
    math::Point p1 = nodes[0].getPoint();
    math::Point p2 = nodes[1].getPoint();


    double pnt1[3] = {p1.X(), p1.Y(), p1.Z()};
    double pnt2[3] = {p2.X(), p2.Y(), p2.Z()};

    double plan0[3] = {APl.getPoint().X(), APl.getPoint().Y(), APl.getPoint().Z()};
    math::Vector3d const &pl_normal = APl.getNormal();
    math::Vector3d pl_e1(pl_normal.Z()-pl_normal.Y(),pl_normal.X()-pl_normal.Z(),pl_normal.Y()-pl_normal.X());
    math::Vector3d pl_e2 = pl_normal.cross(pl_e1);
    math::Point p1_simon = APl.getPoint()+(pl_e1);
    double plan1_simon[3] = {p1_simon.X(),p1_simon.Y(),p1_simon.Z()};
    math::Point p2_simon = APl.getPoint()+(pl_e2);
    double plan2_simon[3] = {p2_simon.X(),p2_simon.Y(),p2_simon.Z()};


    // NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    GEO::Sign sign_1 = GEO::PCK::orient_3d(pnt1, plan0, plan1_simon, plan2_simon);
    GEO::Sign sign_2 = GEO::PCK::orient_3d(pnt2, plan0, plan1_simon, plan2_simon);


    bool intersect_sharp = ((sign_1 > 0 && sign_2 <0) || (sign_1 < 0 && sign_2 >0) );

    if(intersect_sharp){
        //point on the edge
        double Aw0,Aw1;
        int result = APl.intersect(math::Segment(p1,p2), PI,Aw0,Aw1, false) ? 1 : 0;

        //correction
        if(result==0){

            //means we have numerical issues
            // so select p1 or p2
            math::Point proj1 = APl.project(p1);
            math::Point proj2 = APl.project(p2);
            if(p1.distance2(proj1) < p2.distance2(proj2)){
                //P1 in the plane
                //std::cout<<"intersect on p1 "<<std::endl;
                if(onVertexCorrection(nodes[0])) {
                    return -1;
                }

            }
            else{
                //p2 in the plane
                //std::cout<<"intersect on p2 "<<std::endl;
                if(onVertexCorrection(nodes[1])){
                    return -2;
                }
            }
        }
        else{
            double distance1 = PI.distance(p1);
            double distance2 = PI.distance(p2);

            if(distance1 < min_length){
                if(onVertexCorrection(nodes[0])) {
                    //std::cout<<"distance to short with P1"<<std::endl;
                    return -1;
                }
                if(PI.X() == p1.X() && PI.Y() == p1.Y() && PI.Z() == p1.Z()){
                    return 0;
                }
            }

            if(distance2 < min_length){
                //std::cout<<"distance to short with P2"<<std::endl;
                if(onVertexCorrection(nodes[1])) {
                    return -2;
                }
                if(PI.X() == p2.X() && PI.Y() == p2.Y() && PI.Z() == p2.Z()){
                    return 0;
                }
            }

            //We are right well inside the edge
            return 1;
        }
    }
    else if (sign_1==0 && sign_2==0){
        //edge (p1,p2) in the plane
        //std::cout<<"edge on the plane"<<std::endl;
        return 0;
    }
    else if (sign_1==0){
        //p1 in the plane
        //std::cout<<"predicate on p1 "<<std::endl;
        if(onVertexCorrection(nodes[0])) {
            return -1;
        }

    }
    else if (sign_2==0){
        //p2 in the plane
        //std::cout<<"predicate on p2 "<<std::endl;
        if(onVertexCorrection(nodes[1])) {
            return -2;
        }
    }

    return 0;

}
/*----------------------------------------------------------------------------*/

bool DualSheetCreator::onVertexCorrection(Node &ANode) {

    std::vector<Edge> edges = m_mesh->get<Node>(ANode.id()).get<Edge>();

    for (auto const &e : edges) {
        if (m_mesh->isMarked(e, wave_precedent)) {
            //std::cout << "\t --> on wave_prec edge" << std::endl;
            return false;
        }
    }
    return true;
}
/*----------------------------------------------------------------------------*/

const math::Vector3d DualSheetCreator::meanVector(std::vector<gmds::math::Vector3d> AVectors)
{

    math::Vector3d Vfinal = AVectors[0];

    for(int i = 1; i < AVectors.size(); i++){
        Vfinal += AVectors[i];
    }
    Vfinal = Vfinal/AVectors.size();

    Vfinal.normalize();

    return Vfinal;

}
/*----------------------------------------------------------------------------*/

const math::Point DualSheetCreator::meanPoints(std::vector<math::Point> const &APoints){

    return math::Point::massCenter(APoints);
}
/*----------------------------------------------------------------------------*/

double DualSheetCreator::averageInfo(std::vector<intersectInfo> AInfos, intersectInfo &AOutput){

    intersectInfo average;

    std::vector<math::Point> points;
    std::vector<math::Vector3d> vectors;
    std::vector<math::Vector3d> directions;
    for(auto const &info : AInfos){
        points.push_back(info.p);
        vectors.push_back(info.v);
        directions.push_back(info.d);
    }

    average.p = meanPoints(points);
    average.v = meanVector(vectors);
    average.v.normalize();
    average.d = meanVector(directions);
    average.e = AInfos[0].e;
    average.wave = propagation_round;

    AOutput = average;
    double distance = 0;

    for(auto const &info : AInfos){
        distance += average.p.distance(info.p);
    }
    distance = distance / AInfos.size();

    return distance;
}
/*----------------------------------------------------------------------------*/

std::vector<DualSheetCreator::intersectInfo>
DualSheetCreator::face_cut(TCellID ATetra, std::vector<intersectInfo> const &AInfos, bool &move_node, TCellID &moved_id) {

    std::vector<intersectInfo> result;

    std::vector<Face> faces = m_mesh->get<Region>(ATetra).get<Face>();

    for (auto const &info : AInfos) {
        Edge info_edge = m_mesh->get<Edge>(info.e);

        for (auto const &f : faces) {
            std::vector<Edge> edges = f.get<Edge>();
            if (edges[0].id() == info.e || edges[1].id() == info.e || edges[2].id() == info.e) {

                //Test if the face is already treated
                if (!m_mesh->isMarked(f, face_treated)) {


                    //if not, we count the intersection point on the face
                    int nb_inter_points = 0;
                    for (auto const &e : edges) {
                        if (!(*var_info)[e.id()].empty()) {
                            nb_inter_points++;
                        }
                    }

                    //if two or more edges are already intersected, the face is count as treated
                    if (nb_inter_points == 2) {
                        m_mesh->mark(f, face_treated);
                        continue;
                    } else if (nb_inter_points > 2) {
                        continue;
                    }


                    //trying to cut the face
                    intersectInfo new_info = RK_4(f, info, move_node, moved_id);
                    if(new_info.wave != -1){
                        result.push_back(new_info);
                    }
                    if (move_node) {
                        //if the cut has met a node, we return an empty vector
                        result.clear();
                        return result;
                    }

                    m_mesh->mark(f, face_treated);
                }
            }
        }
    }
    return result;
}
/*----------------------------------------------------------------------------*/

DualSheetCreator::intersectInfo
DualSheetCreator::RK_4(const Face &ATriangle, intersectInfo Ainfo, bool &move_node, TCellID &moved_id)
{

    //std::cout<<"================ Runge-kutta ================"<<std::endl;
    intersectInfo result_info;
    result_info.wave = -1;
    std::vector<Edge> edges = ATriangle.get<Edge>();
    std::vector<Node> nodes = ATriangle.get<Node>();

    //init
    math::Point si = Ainfo.p;
    math::Plane A_plane(si, Ainfo.v);

    for(auto e : edges) {
        if (e.id() != Ainfo.e) {

            math::Point B_point;
            int intersect = intersectEdge(A_plane, B_point, e);
            if (intersect > 0) {

                math::Vector3d B = findNormal(B_point, Ainfo.v, e);
                B.normalize();
                math::Plane B_plane(si, B);

                math::Point C_point;
                for(auto eA : edges) {
                    if (eA.id() != Ainfo.e) {
                        intersect = intersectEdge(B_plane, C_point, eA);
                        if (intersect > 0) {
                            math::Vector3d C = findNormal(C_point, Ainfo.v, eA);
                            C.normalize();
                            math::Plane C_plane(si, C);

                            math::Point D_point;
                            for(auto eB : edges) {
                                if (eB.id() != Ainfo.e) {
                                    intersect = intersectEdge(C_plane, D_point, eB);
                                    if (intersect > 0) {

                                        math::Vector3d D = findNormal(D_point, Ainfo.v, eB);
                                        D.normalize();
                                        math::Plane D_plane(si, D);

                                        math::Vector3d vec_result((Ainfo.v + 2 * B + 2 * C + D) / 6);
                                        vec_result.normalize();
                                        math::Plane plane_result(si,vec_result);

                                        math::Point result;
                                        for(auto eC : edges) {
                                            if (eC.id() != Ainfo.e) {
                                                intersect = intersectEdge(plane_result, result, eC);
                                                if (intersect > 0) {

                                                    math::Vector3d propagation_direction(Ainfo.p,result);
                                                    propagation_direction.normalize();
                                                    Ainfo.d.normalize();
                                                    double scalar = propagation_direction.dot(Ainfo.d);

                                                    if(scalar > 0) {

                                                        result_info.p = result;
                                                        result_info.v = findNormal(result, Ainfo.v, eC);
                                                        result_info.v.normalize();
                                                        result_info.tet = Ainfo.tet;
                                                        result_info.wave = propagation_round;
                                                        result_info.e = eC.id();

                                                        //création du vecteur de direction de propagation
                                                        math::Plane plane(result_info.p, result_info.v);

                                                        math::Vector3d vec(Ainfo.p, result_info.p);
                                                        vec.normalize();
                                                        result_info.d = vec;


                                                    }else {
                                                        result_info.wave = -1;
                                                        //exit(40);
                                                    }
                                                    return result_info;

                                                } else if (intersect == -1) {
                                                    move_node = true;
                                                    moved_id = eC.get<Node>()[0].id();
                                                    return result_info;
                                                } else if (intersect == -2) {
                                                    move_node = true;
                                                    moved_id = eC.get<Node>()[1].id();
                                                    return result_info;
                                                }
                                            }
                                        }

                                    } else if (intersect == -1) {
                                        move_node = true;
                                        moved_id = eB.get<Node>()[0].id();
                                        return result_info;
                                    } else if (intersect == -2) {
                                        move_node = true;
                                        moved_id = eB.get<Node>()[1].id();
                                        return result_info;
                                    }
                                }
                            }

                        } else if (intersect == -1) {
                            move_node = true;
                            moved_id = eA.get<Node>()[0].id();
                            return result_info;
                        } else if (intersect == -2) {
                            move_node = true;
                            moved_id = eA.get<Node>()[1].id();
                            return result_info;
                        }
                    }
                }

            } else if (intersect == -1) {
                move_node = true;
                moved_id = e.get<Node>()[0].id();
                return result_info;
            } else if (intersect == -2) {
                move_node = true;
                moved_id = e.get<Node>()[1].id();
                return result_info;
            }
        }
    }
}
/*----------------------------------------------------------------------------*/

bool DualSheetCreator::isAxisSingularity(TCellID AId, const math::Vector3d& AV){

    if((*m_singularities)[AId] > 0) {
        math::Vector3d sing_vec;
        int sing = isSingular(AId, sing_vec);
        // Tetra singular in all axis
        if (sing == 3) {
            std::cout << "singularity in the axis in tet " << AId << std::endl;
            return true;
        }
        else if (sing == 2){

            double cos = AV.dot(sing_vec);

            if (fabs(cos) <= 0.7) {
                std::cout << "singularity direction too far from the direction in tet " << AId << std::endl;
                return true;
            }
        }
    }


    return false;
}
/*----------------------------------------------------------------------------*/

int DualSheetCreator::isSingular(TCellID AId, math::Vector3d& AVOut)
{
    std::vector<TCellID> nodes = m_mesh->get<Region>(AId).getIDs<Node>();
    int sing_type = 0;

    math::Chart charts[4];
    charts[0] = (*m_vertex_chart)[nodes[0]];
    charts[1] = (*m_vertex_chart)[nodes[1]];
    charts[2] = (*m_vertex_chart)[nodes[2]];
    charts[3] = (*m_vertex_chart)[nodes[3]];

    int matchingTab[4][4][3];
    for (unsigned int i = 0; i < 4; i++){
        for (unsigned int j = 0; j < 4; j++){
            charts[i].matchVectors(charts[j],matchingTab[i][j]);
        }
    }

    for (unsigned int i = 0; i < 3; i++){//We take the vectors of charts[0]
        bool goodMatching  = true;
        for (unsigned int j = 1; j < 4; j++){
            for (unsigned int k = 1; k < 4; k++){
                if (matchingTab[j][k][matchingTab[0][j][i]] != matchingTab[0][k][i]){
                    goodMatching = false;
                }
            }
        }
        //
        if (goodMatching){
            AVOut = charts[0][i];
        }
        else{
            sing_type++;
        }
    }

    return sing_type;
}
/*----------------------------------------------------------------------------*/

void DualSheetCreator::clear() {

    for(auto const &t : surface){
        (*m_propagation_round)[t.first] = 1;
    }

    m_mesh->unmarkAll<Edge>(wave_precedent);

    m_mesh->unmarkAll<Face>(face_treated);

    m_mesh->unmarkAll<Region>(mark_wave_tet);
}
/*----------------------------------------------------------------------------*/

void DualSheetCreator::sheet_cleaning() {

    clear();

    //for(auto f : m_mesh->faces()){

    //if(f == 2201) exit(0);
    //}
    //std::cout<<2201<<":"<<linker.getGeomDim<Face>(2201)<<std::endl;

    int nb_edge = 2;

    std::cout<<"=============== Cleaning sheet =============== "<<std::endl;

    std::map<gmds::TCellID,gmds::math::Vector3d> surface_it = surface;
    std::map<gmds::TCellID,gmds::math::Vector3d> surface_tmp;

    bool modification = false;

    do {

        m_mesh->unmarkAll<Edge>(wave_precedent);
        m_mesh->unmarkAll<Region>(mark_wave_tet);

        modification = false;
        for (auto const &t : surface_it) {
            std::vector<Edge> edges = m_mesh->get<Region>(t.first).get<Edge>();
            for (auto const &e : edges) {

                int changement = 0;
                int color;
                if (!m_mesh->isMarked(e, wave_precedent)) {
                    std::vector<Region> tets = e.get<Region>();
                    Region tet;
                    if (linker.getGeomDim(e) == 3) {
                        bool breaker = false;
                        for (auto const &tet_bounder : tets) {
                            std::vector<TCellID> faces = tet_bounder.getIDs<Face>();
                            for (auto const &f : faces) {
                                if (linker.getGeomDim<Face>(f) == 3) {
                                    //std::cout<<"tet on surf"<<std::endl;
                                    tet = tet_bounder;
                                    breaker = true;
                                }
                            }
                            if (breaker) break;
                        }
                    } else {
                        tet = tets[0];
                    }
                    color = (*m_propagation_round)[tet.id()];


                    while (!m_mesh->isMarked(tet, mark_wave_tet)) {

                        m_mesh->mark<Region>(tet.id(), mark_wave_tet);

                        std::vector<Face> faces = tet.get<Face>();
                        for (auto const &f : faces) {
                            std::vector<Region> tets_of_face = f.get<Region>();
                            if (!m_mesh->isMarked(tets_of_face[0], mark_wave_tet)) {

                                if (std::find(tets.begin(), tets.end(), tets_of_face[0]) != tets.end()) {
                                    tet = tets_of_face[0];
                                    if (color != (*m_propagation_round)[tet.id()]) {
                                        color = (*m_propagation_round)[tet.id()];
                                        changement++;
                                    }

                                    break;
                                }
                            } else if (tets_of_face.size() > 1) {
                                if (!m_mesh->isMarked(tets_of_face[1], mark_wave_tet)) {
                                    if (std::find(tets.begin(), tets.end(), tets_of_face[1]) != tets.end()) {
                                        tet = tets_of_face[1];
                                        if (color != (*m_propagation_round)[tet.id()]) {
                                            color = (*m_propagation_round)[tet.id()];
                                            changement++;
                                        }

                                        break;
                                    }
                                }
                            }
                        }
                    }


                    if ((changement > 2) || (changement > 1 && linker.getGeomDim(e) == 3)) {

                        modification = true;

                        //std::cout<<"Edge non valid"<<std::endl;

                        math::Vector3d vec;
                        for (auto const &te : tets) {
                            if ((*m_propagation_round)[te.id()] == 1) {
                                vec = surface[te.id()];
                            }
                        }
                        for (auto const &te : tets) {
                            if ((*m_propagation_round)[te.id()] == -1) {
                                //std::cout<<"tet added "<<te.id()<<std::endl;

                                surface.emplace(te.id(), vec);
                                (*m_propagation_round)[te.id()] = 1;
                                surface_tmp.emplace(te.id(), vec);

                            }
                        }
                        nb_edge++;
                    }

                    for (auto const &tet_clean : tets) {
                        m_mesh->unmark(tet_clean, mark_wave_tet);
                    }

                    m_mesh->mark(e, wave_precedent);
                }
            }
        }
        surface_it = surface_tmp;
    }while(modification);

    std::cout<<"End cleaning"<<std::endl;
}
/*----------------------------------------------------------------------------*/

void DualSheetCreator::reset() {
    clear();
    for(auto const &t : surface){
        (*m_propagation_round)[t.first] = 0;
        std::vector<TCellID> edges = m_mesh->get<Region>(t.first).getIDs<Edge>();
        for(auto e : edges){
            (*var_info)[e].clear();
        }
    }

}
/*----------------------------------------------------------------------------*/




























