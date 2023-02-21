/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/smoothy/AngleBasedQuadSmoother.h>
#include <gmds/smoothy/LaplacianSmoother.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/blockMesher/BlockMesher.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== Blocker ================" << std::endl;

//==================================================================
// PARAMETERS' PARSING
//==================================================================
    std::string file_geom, file_mesh, file_mesh_out, file_block_out;
    int nb_curve_smooth_iterations=0;
    int nb_surface_smooth_iterations=0;
    int nb_volume_smooth_iterations=0;
    int edge_discretization=0;
    if (argc != 9) {
        std::cout << "Require four paramaters : \n";
        std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
        std::cout << "  - [IN ] mesh (.vtk) that describes the block structure to refine and smoothed, \n";
        std::cout << "  - [IN ] edge discretization of each block edge,\n";
        std::cout << "  - [IN ] number of curve smoothing iterations,\n";
        std::cout << "  - [IN ] number of surface smoothing iterations,\n";
        std::cout << "  - [IN ] number of volume smoothing iterations,\n";
        std::cout << "  - [OUT] the final blocks (.vtk),\n";
        std::cout << "  - [OUT] the final mesh (.vtk). \n"<< std::endl;
        throw gmds::GMDSException("Wrong number of parameters");
    }

    file_geom = std::string(argv[1]);
    file_mesh = std::string(argv[2]);
    nb_curve_smooth_iterations = atoi(std::string(argv[3]).c_str());
    nb_surface_smooth_iterations = atoi(std::string(argv[4]).c_str());
    nb_volume_smooth_iterations = atoi(std::string(argv[5]).c_str());
    edge_discretization = atoi(std::string(argv[6]).c_str());
    file_block_out = std::string(argv[7]);
    file_mesh_out = std::string(argv[8]);
    std::cout << "Parameters " << std::endl;
    std::cout << "  - Geometry file: " << file_geom << std::endl;
    std::cout << "  - Mesh file    : " << file_mesh << std::endl;
    std::cout << "  - Nb curve iterations: " << nb_curve_smooth_iterations << std::endl;
    std::cout << "  - Nb surface iterations: " << nb_surface_smooth_iterations << std::endl;
    std::cout << "  - Nb volume iterations: " << nb_volume_smooth_iterations << std::endl;
    std::cout << "  - Edge discretization: " << edge_discretization << std::endl;
    std::cout << "  - Output block file  : " << file_block_out << std::endl;
    std::cout << "  - Output mesh file  : " << file_mesh_out << std::endl;
    std::cout << "=======================================" << std::endl;

//==================================================================
// GEOMETRY READING
//==================================================================
    std::cout<<"> Start geometry reading"<<std::endl;
    Mesh geometry(MeshModel(DIM3 | R | F | E | N |
                            R2N | R2F | R2E |
                            F2N | F2R | F2E |
                            E2F | E2N | N2E | N2R | N2F));

    IGMeshIOService ioService(&geometry);
    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(N| R);
    vtkReader.read(file_geom);
    MeshDoctor doc(&geometry);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager geom_manager;
    geom_manager.initFrom3DMesh(&geometry);

//==================================================================
// MESH READING
//==================================================================
    std::cout<<"> Start block reading"<<std::endl;
//the used model is specified according to the geom smoother requirements.
    Mesh blocking(MeshModel(DIM3 | R | F | E | N |
                            R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

    IGMeshIOService ioService2(&blocking);
    VTKReader vtkReader2(&ioService2);
    vtkReader2.setCellOptions(N|R);
    vtkReader2.read(file_mesh);
    MeshDoctor doc2(&blocking);
    doc2.buildFacesAndR2F();
    doc2.buildEdgesAndX2E();
    doc2.updateUpwardConnectivity();

//==================================================================
// MARK ALL THE BOUNDARY CELL OF THE INIT MESH
//==================================================================

    std::cout<<"> Start block boundary retrieval"<<std::endl;
//we get all the nodes that are on the mesh boundary
    BoundaryOperator op(&blocking);
    auto mark_node_NAN = blocking.newMark<Node>();
    auto mark_node_on_pnt = blocking.newMark<Node>();
    auto mark_node_on_crv = blocking.newMark<Node>();
    auto mark_node_on_srf = blocking.newMark<Node>();
    auto mark_edge_on_crv = blocking.newMark<Edge>();
    auto mark_edge_on_srf = blocking.newMark<Edge>();
    auto mark_face_on_srf = blocking.newMark<Face>();

    op.markCellOnGeometry(mark_face_on_srf,
                          mark_edge_on_srf,
                          mark_node_on_srf,
                          mark_edge_on_crv,
                          mark_node_on_crv,
                          mark_node_on_pnt,
                          mark_node_NAN);

    IGMeshIOService ioService_geom(&geometry);
    VTKWriter writer_geom(&ioService_geom);
    writer_geom.setCellOptions(N|F);
    writer_geom.setDataOptions(N|F);
    writer_geom.write("surf_input_classification.vtk");
//==================================================================
// TRY TO CLASSIFY MESH CELLS
// The classification process is not straightforward as we have to
// consider that the mesh to classify is a "coarse" version of the
// geometry. For instance, if you try and project the middle point
// of an edge onto a curve, the distance between the initial point
// and the projected one will be probably not equal to zero. While,
// the same point will be lying on adjacent surface with a null
// distance. So the strategy is to
// 1. Classify faces on surfaces
// 2. Get edges that are shared by two faces classified onto 2
//    surfaces. Those edges must be classified on curves. It is
//    also the case of the edge endpoint
// 3. Do the same for nodes shared by at least 3 faces classified
//    onto 3 surfaces. Those nodes must be classified onto points.
//==================================================================

    std::cout<<"> Start block->geometry classification"<<std::endl;
    cad::GeomMeshLinker linker(&blocking, &geom_manager);
    std::vector<cad::GeomSurface*> surfaces;
    std::vector<cad::GeomCurve*> curves;
    std::vector<cad::GeomPoint*> points;

    geom_manager.getSurfaces(surfaces);
    geom_manager.getCurves(curves);
    geom_manager.getPoints(points);

//==================================================================flatpak install flathub md.obsidian.Obsidian
//First, we classify each face

    for(auto f_id: blocking.faces()){
        Face f= blocking.get<Face>(f_id);
        if(blocking.isMarked(f, mark_face_on_srf) ){
            //we've got a boundary face
            math::Point p = f.center();
            math::Vector3d v = f.normal();
            //we ensure to get the outward normal that goes from the inside block
            //towards the outside
            //as f is a boundary face, it is connected to only one block
            Region r = f.get<Region>()[0];
            math::Vector3d v2inside= r.center()-p;
            if(v.dot(v2inside)>0){
                //v points inside
                v = v.opp();
            }

            v.normalize();
            std::cout<<"FACE TO PROJECT "<<f_id<<" with center point  "<<p<<" and normal direction "<<v<<std::endl;
            double min_dist = 100000;
            int min_entity_dim=-1;
            int min_entity_id = -1;
            std::map<TCellID , double> surf_dist;
            std::map<TCellID , double> surf_dot;
            bool found_surf=false;
            for(auto s:surfaces){
                math::Point closest_pnt=p;
                s->project(closest_pnt);
                double dist = p.distance2(closest_pnt);
                math::Vector v_proj=closest_pnt-p;
                v_proj.normalize();
                surf_dot[s->id()]=(v.dot(v_proj));
                surf_dist[s->id()]=dist;

                if(dist<min_dist) {
                    //either the points are the same or we check the vector orientation
                    if (dist<1e-4 || fabs(v.dot(v_proj)) > 0.2 ) {
                        bool keep_projection = true;

                       if (v.dot(v_proj)<0){
                            //we have to check that we do not cross the volume. It can occur for thin curved
                            //shape (like cylinder with a hole)

                            //We take other points closed to the face corners and we project all of them, if they go to
                            // the same surf, we keep it. Otherwise, we will put another surface later.
                            keep_projection=false;
                            std::vector<Node> corners = f.get<Node>();
                            std::vector<math::Point> corner_pnts;
                            for(auto c:corners){
                                corner_pnts.push_back(0.1*p+0.9* c.point());
                            }
                            std::set<int> corner_surf;
                            for(auto cp:corner_pnts) {
                                double cp_min_dist = 100000;
                                auto cp_min_surf_id = -1;
                                for (auto s: surfaces) {
                                    math::Point closest_pnt = cp;
                                    s->project(closest_pnt);

                                    double local_dist = cp.distance2(closest_pnt);
                                    if (local_dist < cp_min_dist) {
                                        cp_min_dist = local_dist;
                                        cp_min_surf_id = s->id();
                                    }
                                }
                                corner_surf.insert(cp_min_surf_id);
                            }
                            if(corner_surf.size()==1 && (*corner_surf.begin()==s->id())){
                                //ok we project on it
                                keep_projection= true;
                            }

                        }
                        if(keep_projection) {
                            min_dist = dist;
                            min_entity_dim = 2;
                            min_entity_id = s->id();
                            found_surf = true;
                        }
                    }
                }
            }
            if(!found_surf){
                min_dist = 100000;
                min_entity_dim=-1;
                min_entity_id = -1;
                //means we have a block face in a concave area
                for(auto s:surfaces){
                    math::Point closest_pnt=p;
                    s->project(closest_pnt);

                    double dist = p.distance2(closest_pnt);
                    math::Vector v_proj=closest_pnt-p;
                    v_proj.normalize();
                    surf_dot[s->id()]=(v.dot(v_proj));
                    surf_dist[s->id()]=dist;
                    if(dist<min_dist) {
                        min_dist = dist;
                        min_entity_dim = 2;
                        min_entity_id = s->id();
                        found_surf=true;

                    }
                }
            }
            if (min_entity_dim==2){
                std::cout<<"==> Link face "<<f_id<<" on surf "<<min_entity_id<<" with distance:" <<min_dist<<std::endl;
                linker.linkFaceToSurface(f_id,min_entity_id);
            }
            else{
                std::cout<<"\t ====> Link error for classifying a face"<<std::endl;
                throw GMDSException("Link error for classifying a face");
            }

        }
    }
//==================================================================
//Second, we classify each edge
    for(auto e_id: blocking.edges()){
        Edge e= blocking.get<Edge>(e_id);
        if(blocking.isMarked(e, mark_edge_on_crv) ||
           blocking.isMarked(e, mark_edge_on_srf) ){
//we've got a boundary edge, now we get the 2 boundary faces
// around e
            std::vector<Node> e_nodes = e.get<Node>();
            std::vector<TCellID> adj_face_id = e.getIDs<Face>();
            std::vector<TCellID> adj_bnd_face_ids;
            for(auto f_id:adj_face_id){
                if(blocking.isMarked<Face>(f_id, mark_face_on_srf)){
                    adj_bnd_face_ids.push_back(f_id);
                }
            }
            if(adj_bnd_face_ids.size()!=2){
                std::cout<<"ERROR: One boundary edge is adjacent to more "
                         <<"than 2 boundary faces ("
                         <<"face "<<e_id<<")"<<std::endl;
            }
            auto surf0_id = linker.getGeomId<Face>(adj_bnd_face_ids[0]);
            auto surf1_id = linker.getGeomId<Face>(adj_bnd_face_ids[1]);

            if(surf0_id==surf1_id){
//edge embedded in the surface
                linker.linkEdgeToSurface(e_id,surf1_id);
            }
            else{
                //it must be an edge classified on curves
                // Warning, like for the face, it is not necessary the geometrically closest to be connected to.
                //We check which curve is bounded by surf0_id and surf1_id
                bool found = false;

                std::vector<cad::GeomCurve*> candidate_curves;
                for(auto ci:curves) {
                    std::vector<cad::GeomSurface *> c_surfs = ci->surfaces();
                    if (c_surfs.size() == 2) {
                        if ((c_surfs[0]->id() == surf0_id && c_surfs[1]->id() == surf1_id) ||
                            (c_surfs[1]->id() == surf0_id && c_surfs[0]->id() == surf1_id)) {
                            std::cout<<"Curve "<<ci->id()<<" adj to surfaces "<<c_surfs[0]->id()<<" and "<<c_surfs[1]->id()<<std::endl;
                            found = true;
                            candidate_curves.push_back(ci);
                        }
                    }
                }
                if(found) {
                    if(candidate_curves.size()==1) {
                        linker.linkEdgeToCurve(e_id, candidate_curves[0]->id());
                    }
                    else{
                        //we project on each of curve and keep the closest one
                        double min_dist = 1e6;
                        cad::GeomCurve* selected_curve = NULL;
                        math::Point center_edge = e.center();
                        for(auto ci:candidate_curves) {
                            math::Point pi = center_edge;
                            ci->project(pi);
                            auto di = pi.distance2(center_edge);
                            if(di<min_dist){
                                min_dist=di;
                                selected_curve = ci;
                            }
                        }

                        linker.linkEdgeToCurve(e_id, selected_curve->id());

                    }
                }
                else{
                    throw GMDSException("BlockMesher error: impossible to link a block edge onto the geometry");
                }
            }
        }
    }
//==================================================================
//we classify each node
    auto on_pnt=0, on_curve=0, on_surf=0;
    for(auto n_id: blocking.nodes()){
        Node n= blocking.get<Node>(n_id);
        if(blocking.isMarked(n, mark_node_on_pnt) ||
           blocking.isMarked(n, mark_node_on_crv) ||
           blocking.isMarked(n, mark_node_on_srf) ){
//we've got a boundary node
            //As face and edges are already classified, we used it
            // A node is on a surface if all its adjacent boundary faces are on the same
            // surface. If they are on two it is on a curve potentially
            std::vector<TCellID> adj_faces = n.getIDs<Face>();
            std::vector<TCellID> adj_bnd_faces;
            for(auto id_face:adj_faces){
                if(blocking.isMarked<Face>(id_face, mark_face_on_srf)) {
                    adj_bnd_faces.push_back(id_face);
                }
            }
            std::set<int> bnd_surfaces;

            for (auto  id_face:adj_bnd_faces){
                bnd_surfaces.insert(linker.getGeomId<Face>(id_face));
            }
            int min_entity_dim=-1;
            int min_entity_id = -1;
            double min_dist = 100000;

            if(bnd_surfaces.size()==1){
                //on a single surface
                min_entity_dim=2;
                min_entity_id=*(bnd_surfaces.begin());
            }
            else if(bnd_surfaces.size()==2){
                //on a curve or a point that is connected to a single curve
                math::Point node_loc = n.point();
                for(auto p:points){
                    math::Point closest_pnt=node_loc;
                    double dist = node_loc.distance(p->point());
                    if(dist<min_dist){
                        min_dist =dist;
                        min_entity_dim=0;
                        min_entity_id=p->id();
                    }
                }
                for(auto c:curves){
                    math::Point closest_pnt=node_loc;
                    c->project(closest_pnt);
                    double dist = node_loc.distance(closest_pnt);
                    if(dist<min_dist){
//WARNING: Take care of this trick that is not good at all but mandatory to be staying on the
// curve and not on the surface
                        if(dist<1e-4)
                            dist=0;
                        min_dist =dist;
                        min_entity_dim=1;
                        min_entity_id=c->id();
                    }
                }

            } else{
                //On a point!!
                math::Point node_loc = n.point();
                for(auto p:points){
                    math::Point closest_pnt=node_loc;
                    double dist = node_loc.distance(p->point());
                    if(dist<min_dist){
                        min_dist =dist;
                        min_entity_dim=0;
                        min_entity_id=p->id();
                    }
                }
            }

            if(min_entity_dim==0){
                on_pnt++;
                linker.linkNodeToPoint(n_id,min_entity_id);
                std::cout<<"Node "<<n_id<<" is on point "<<min_entity_id<<std::endl;
            }
            else if (min_entity_dim==1){
                on_curve++;
                linker.linkNodeToCurve(n_id,min_entity_id);
            }
            else if (min_entity_dim==2){
                on_surf++;
                linker.linkNodeToSurface(n_id,min_entity_id);
            }
            else{
                throw GMDSException("Link error for classifying a node");
            }

        }
    }
    std::cout<<"  info [node classified on points "<<on_pnt;
    std::cout<<", on curves "<<on_curve;
    std::cout<<", on surfs "<<on_surf<<"]"<<std::endl;

    linker.writeVTKDebugMesh("linker_debug.vtk");
//==================================================================
// PERFORM THE BLOCK SMOOTHING NOW
//==================================================================
    std::cout<<"> Start block smoothing"<<std::endl;
    smoothy::LaplacianSmoother smoother(&linker);
    if(!smoother.isValid())
    {
        std::cout<<"INVALID MODEL"<<std::endl;
        exit(1);
    }
    std::cout<<"  - start  smoothing ("<<nb_curve_smooth_iterations<<" iterations)"<<std::endl;
    smoother.smoothCurves(nb_curve_smooth_iterations);
    smoother.smoothSurfaces(nb_surface_smooth_iterations);
   // std::cout<<"  - start volume smoothing ("<<nb_volume_smooth_iterations<<" iterations)"<<std::endl;
   // smoother.smoothVolumes(nb_volume_smooth_iterations);

//==================================================================
// COMPUTE BLOCK QUALITY
//==================================================================
    Variable<double>* var_quality = blocking.newVariable<double,gmds::GMDS_REGION>("GMDS_Scaled_Jacobian");
    for(auto r_id:blocking.regions()){
        Region r = blocking.get<Region>(r_id);
        var_quality->value(r_id)=fabs(r.computeScaledJacobian());
    }
    std::cout << "> Write smoothed block mesh in: " << file_mesh_out << std::endl;
    IGMeshIOService ioService3(&blocking);
    VTKWriter vtkWriter(&ioService3);
    vtkWriter.setCellOptions(N|R);
    vtkWriter.setDataOptions(N|R);
    vtkWriter.write(file_block_out);

//==================================================================
// NOW WE CREATE THE FINAL MESH
//==================================================================
    //The first stage consists in performing a transfinite interpolation of blocks
    BlockMesher bm(&blocking,&linker);
    std::cout<<"> Start mesh generation via transfinite interpolation"<<std::endl;
    bm.execute(edge_discretization);
    //All boundary node, edge and curves are classified.
    std::map<TCellID, Cell::Data> mesh_node_info = bm.mesh_node_classification();
    std::map<TCellID, Cell::Data> mesh_edge_info = bm.mesh_edge_classification();
    std::map<TCellID, Cell::Data> mesh_face_info = bm.mesh_face_classification();

    std::cout<<"> Start block smoothing"<<std::endl;
    bm.mesh()->changeModel(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|N2F|N2E|N2R));
    MeshDoctor doc3(bm.mesh());
    doc3.updateUpwardConnectivity();
    cad::GeomMeshLinker linker_final(bm.mesh(),&geom_manager);
    for(auto info:mesh_node_info){
        Node n = bm.mesh()->get<Node>(info.first);
        if(info.second.dim==0) {
            linker_final.linkNodeToPoint(n.id(), info.second.id);
        }
        else if(info.second.dim==1) {
            linker_final.linkNodeToCurve(n.id(), info.second.id);
        }
        else if(info.second.dim==2)
            linker_final.linkNodeToSurface(n.id(),info.second.id);
    }
    for(auto info:mesh_edge_info){
        Edge e = bm.mesh()->get<Edge>(info.first);
        if(info.second.dim==1) {
            linker_final.linkEdgeToCurve(e.id(), info.second.id);
            std::cout<<"Edge "<<e.getIDs<Node>()[0]<<"-"<<e.getIDs<Node>()[1]<<" on curve "<<info.second.id<<std::endl;
        }
    }
    for(auto info:mesh_face_info){
        Face f = bm.mesh()->get<Face>(info.first);
        if(info.second.dim==2)
            linker_final.linkFaceToSurface(f.id(),info.second.id);
    }
    smoothy::LaplacianSmoother smoother_final(&linker_final);
    if(!smoother_final.isValid())
    {
        std::cout<<"INVALID MODEL"<<std::endl;
        exit(1);
    }
    std::cout<<"  - start smoothing ("<<nb_curve_smooth_iterations<<" iterations)"<<std::endl;
    smoother_final.smoothCurves(nb_curve_smooth_iterations);
    smoother_final.smoothSurfaces(nb_surface_smooth_iterations);
    std::cout<<"> Start final mesh writing"<<std::endl;
    IGMeshIOService ioService_mesh(bm.mesh());
    VTKWriter vtkWriter_mesh(&ioService_mesh);
    vtkWriter_mesh.setCellOptions(N|F);
    vtkWriter_mesh.setDataOptions(N|F);
    vtkWriter_mesh.write(file_mesh_out);

    VTKWriter vtkWriter_curves(&ioService_mesh);
    vtkWriter_curves.setCellOptions(N|E);
    vtkWriter_curves.setDataOptions(N|E);
    vtkWriter_curves.write("mesh_curves.vtk");

    std::cout << "======== Task done by blocker =========" << std::endl;
}