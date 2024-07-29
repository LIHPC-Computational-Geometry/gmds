/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/smoothy/LaplacianSmoother3C.h>
/*----------------------------------------------------------------------------*/
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "============== Smoothy ================" << std::endl;

    //==================================================================
    // PARAMETERS' PARSING
    //==================================================================
    std::string file_geom, file_mesh, file_out;
    int nb_iterations=0;
    if (argc != 5) {
        std::cout << "Require four paramaters : \n";
        std::cout << "  - [IN ] tetrahedral mesh (.vtk) that describes the geometry, \n";
        std::cout << "  - [IN ] mesh (.vtk) that describes the mesh to be smoothed, \n";
        std::cout << "  - [IN ] number of smoothing iterations,\n";
        std::cout << "  - [OUT] the smoothed mesh (.vtk). \n"<< std::endl;
        throw gmds::GMDSException("Wrong number of parameters");
    }

    file_geom = std::string(argv[1]);
    file_mesh = std::string(argv[2]);
    nb_iterations = atoi(std::string(argv[3]).c_str());
    file_out = std::string(argv[4]);
    std::cout << "Parameters " << std::endl;
    std::cout << "  - Geometry file: " << file_geom << std::endl;
    std::cout << "  - Mesh file    : " << file_mesh << std::endl;
    std::cout << "  - Nb iterations: " << nb_iterations << std::endl;
    std::cout << "  - Output file  : " << file_out << std::endl;
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
    std::cout<<"> Start mesh reading"<<std::endl;
    //the used model is specified according to the geom smoother requirements.
    Mesh m(MeshModel(DIM3 | R | F | E | N |
                     R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | N2F));

    IGMeshIOService ioService2(&m);
    VTKReader vtkReader2(&ioService2);
    vtkReader2.setCellOptions(N|R);
    vtkReader2.read(file_mesh);
    MeshDoctor doc2(&m);
    doc2.buildFacesAndR2F();
    doc2.buildEdgesAndX2E();
    doc2.updateUpwardConnectivity();

    //==================================================================
    // MARK ALL THE BOUNDARY CELL OF THE INIT MESH
    //==================================================================

    std::cout<<"> Start mesh boundary retrieval"<<std::endl;
    //we get all the nodes that are on the mesh boundary
    BoundaryOperator op(&m);
    auto mark_node_NAN = m.newMark<Node>();
    auto mark_node_on_pnt = m.newMark<Node>();
    auto mark_node_on_crv = m.newMark<Node>();
    auto mark_node_on_srf = m.newMark<Node>();
    auto mark_edge_on_crv = m.newMark<Edge>();
    auto mark_edge_on_srf = m.newMark<Edge>();
    auto mark_face_on_srf = m.newMark<Face>();

    op.markCellOnGeometry(mark_face_on_srf,
                          mark_edge_on_srf,
                          mark_node_on_srf,
                          mark_edge_on_crv,
                          mark_node_on_crv,
                          mark_node_on_pnt,
                          mark_node_NAN);


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

    std::cout<<"> Start mesh->geometry classification"<<std::endl;
    cad::GeomMeshLinker linker(&m, &geom_manager);
    std::vector<cad::GeomSurface*> surfaces;
    std::vector<cad::GeomCurve*> curves;
    std::vector<cad::GeomPoint*> points;

    geom_manager.getSurfaces(surfaces);
    geom_manager.getCurves(curves);
    geom_manager.getPoints(points);

    //==================================================================
    //First, we classify each face
    for(auto f_id: m.faces()){
        Face f= m.get<Face>(f_id);
        if(m.isMarked(f, mark_face_on_srf) ){
            //we've got a boundary face
            math::Point p = f.center();

            double min_dist = 100000;
            int min_entity_dim=-1;
            int min_entity_id = -1;
            for(auto s:surfaces){
                math::Point closest_pnt=p;
                s->project(closest_pnt);
                double dist = p.distance2(closest_pnt);
                if(dist<min_dist){
                    min_dist =dist;
                    min_entity_dim=2;
                    min_entity_id=s->id();
                }
            }
            if (min_entity_dim==2){
                linker.linkFaceToSurface(f_id,min_entity_id);
            }
            else{
                throw GMDSException("Link error for classifying a face");
            }
        }
    }
    //==================================================================
    //Second, we classify each edge
    for(auto e_id: m.edges()){
        Edge e= m.get<Edge>(e_id);
        if(m.isMarked(e, mark_edge_on_crv) ||
        m.isMarked(e, mark_edge_on_srf) ){
            //we've got a boundary edge, now we get the 2 boundary faces
            // around e
            std::vector<Node> e_nodes = e.get<Node>();
            std::vector<TCellID> adj_face_id = e.getIDs<Face>();
            std::vector<TCellID> adj_bnd_face_ids;
            for(auto f_id:adj_face_id){
                if(m.isMarked<Face>(f_id,mark_face_on_srf)){
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
                //we look for the closest curve now
                math::Point p0 = e_nodes[0].point();
                math::Point p1 = e_nodes[1].point();
                math::Point p = 0.5*(p0+p1);

                double min_dist = 100000;
                int min_entity_dim=-1;
                int min_entity_id = -1;
                for(auto c:curves){
                    math::Point closest_pnt=p;
                    c->project(closest_pnt);
                    double dist = p.distance2(closest_pnt);
                    if(dist<min_dist){
                        min_dist =dist;
                        min_entity_dim=1;
                        min_entity_id=c->id();
                    }
                }
                linker.linkEdgeToCurve(e_id,min_entity_id);
                //std::cout<<"Edge "<<e_nodes[0].id()<<"-"<<e_nodes[1].id()<<" is on curve "<<min_entity_id<<std::endl;
            }
        }
    }
    //==================================================================
    //we classify each node
    auto on_pnt=0, on_curve=0, on_surf=0;
    for(auto n_id: m.nodes()){
        Node n= m.get<Node>(n_id);
        if(m.isMarked(n, mark_node_on_pnt) ||
           m.isMarked(n, mark_node_on_crv) ||
           m.isMarked(n, mark_node_on_srf) ){
            //we've got a boundary node
            math::Point node_loc = n.point();
            double min_dist = 100000;
            int min_entity_dim=-1;
            int min_entity_id = -1;
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
            for(auto s:surfaces){
                math::Point closest_pnt=node_loc;
                s->project(closest_pnt);
                double dist = node_loc.distance(closest_pnt);
                if(dist<min_dist){
                    min_dist =dist;
                    min_entity_dim=2;
                    min_entity_id=s->id();
                }
            }

            if(min_entity_dim==0){
                on_pnt++;
                linker.linkNodeToPoint(n_id,min_entity_id);
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
    // PERFORM THE MESH SMOOTHING NOW
    //==================================================================
    std::cout<<"> Start smoothing"<<std::endl;
    smoothy::LaplacianSmoother3C smoother(&m,&linker);
    if(!smoother.isValid())
    {
        std::cout<<"INVALID MODEL"<<std::endl;
        exit(1);
    }
	 smoother.setNbIterations(nb_iterations);
    std::cout<<"  - start curve smoothing ("<<nb_iterations<<" iterations)"<<std::endl;
    smoother.smoothCurves();
    std::cout<<"  - start surface smoothing ("<<nb_iterations<<" iterations)"<<std::endl;
    smoother.smoothSurfaces();
    std::cout<<"  - start volume smoothing ("<<nb_iterations<<" iterations)"<<std::endl;
    smoother.smoothVolumes();

    //==================================================================
    // COMPUTE FINAL QUALITY
    //==================================================================
    Variable<double>* var_quality = m.newVariable<double,gmds::GMDS_REGION>("GMDS_Scaled_Jacobian");
    for(auto r_id:m.regions()){
        Region r = m.get<Region>(r_id);
        var_quality->value(r_id)=fabs(r.computeScaledJacobian());
    }
    std::cout<<"> Write output mesh in: "<<file_out<<std::endl;
    IGMeshIOService ioService3(linker.mesh());
    VTKWriter vtkWriter(&ioService3);
    vtkWriter.setCellOptions(N|R);
    vtkWriter.setDataOptions(N|R);
    vtkWriter.write(file_out);
    std::cout << "======== Task done by Smoothy =========" << std::endl;
}