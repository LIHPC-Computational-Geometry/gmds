/*----------------------------------------------------------------------------*/
/*
 * FACManager.t.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/ig/MeshDoctor.h>

#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <set>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/igalgo/BoundaryExtractor2D.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        FACManager::FACManager()
                :m_mesh(DIM3|N|E|F|F2N|N2F|E2N|E2F|F2E|N2E)
        {
            FACPoint::resetIdCounter();
            FACCurve::resetIdCounter();
            FACSurface::resetIdCounter();
            FACVolume::resetIdCounter();
        }

/*----------------------------------------------------------------------------*/
        FACManager::~FACManager()
        {
            for(unsigned int i=0; i < m_volumes.size(); i++)
                if(m_volumes[i] != nullptr)
                    delete m_volumes[i];

            for(unsigned int i=0; i < m_surfaces.size(); i++)
                if(m_surfaces[i] != nullptr)
                    delete m_surfaces[i];

            for(unsigned int i=0; i < m_curves.size(); i++)
                if(m_curves[i] != nullptr)
                    delete m_curves[i];

            for(unsigned int i=0; i < m_points.size(); i++)
                if(m_points[i] != nullptr)
                    delete m_points[i];

        }
/*----------------------------------------------------------------------------*/
        void FACManager::initFrom3DMesh(Mesh* AFromMesh){


            BoundaryExtractor3D boundary_extractor(AFromMesh, &m_mesh);


            Variable<int>* node_on_pnt = nullptr;

            Variable<int>* node_on_crv = nullptr;
            Variable<int>* node_on_srf = nullptr;
            Variable<int>* edge_on_crv = nullptr;
            Variable<int>* face_on_srf = nullptr;
            try {
                node_on_pnt= m_mesh.newVariable<int, GMDS_NODE>("on_point"  );
                node_on_crv = m_mesh.newVariable<int, GMDS_NODE>("on_curve"  );
                node_on_srf = m_mesh.newVariable<int, GMDS_NODE>("on_surface");
                edge_on_crv = m_mesh.newVariable<int, GMDS_EDGE>("on_curve"  );
                face_on_srf = m_mesh.newVariable<int, GMDS_FACE>("on_surface");

            }
            catch (GMDSException &e){
                node_on_pnt= m_mesh.getVariable<int, GMDS_NODE>("on_point"  );
                node_on_crv = m_mesh.getVariable<int, GMDS_NODE>("on_curve"  );
                node_on_srf = m_mesh.getVariable<int, GMDS_NODE>("on_surface");
                edge_on_crv = m_mesh.getVariable<int, GMDS_EDGE>("on_curve"  );
                face_on_srf = m_mesh.getVariable<int, GMDS_FACE>("on_surface");

            }

            boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, face_on_srf, node_on_srf,node_on_crv);

            boundary_extractor.execute();

            //to update the cad model from m_mesh, we have to add cell groups into m_mesh.

            gmds::MeshDoctor doc(&m_mesh);
            doc.buildEdgesAndX2E();
            doc.updateUpwardConnectivity();
            // vertices
            for(auto n_id: m_mesh.nodes()){
                if(node_on_pnt->value(n_id)!=0){
                    FACPoint* p = new FACPoint(&m_mesh, n_id, node_on_pnt->getName());
                    m_points.push_back(p);
                    m_map_node_var_to_point[node_on_pnt->value(n_id)]=p->id();
                }
            }

            // curves
            std::map<int,std::vector<TCellID> > edges_per_curve;
            std::map<int,std::string> curve_names;
            for(auto e_id: m_mesh.edges()){
                if (edge_on_crv->value(e_id) != 0) {
                    edges_per_curve[edge_on_crv->value(e_id)].push_back(e_id);
                    curve_names[edge_on_crv->value(e_id)] = edge_on_crv->getName();
                }
            }

            for(auto e2c: edges_per_curve) {
                std::vector<TCellID> curve_edges = e2c.second;

                //we get all the nodes of the edge
                std::set<TCellID> all_ids;
                for (auto e_id: curve_edges){
                    Edge e = m_mesh.get<Edge>(e_id);
                    std::vector<TCellID> n_ids = e.getIDs<Node>();
                    all_ids.insert(n_ids[0]);
                    all_ids.insert(n_ids[1]);
                }
                std::vector<TCellID> node_ids;
                node_ids.assign(all_ids.begin(),all_ids.end());
                FACCurve* c =
                        new FACCurve(&m_mesh,node_ids,curve_edges, curve_names[e2c.first]);
                m_curves.push_back(c);

                m_map_edge_var_to_curve[e2c.first]=c->id();
            }

            std::map<int,std::vector<TCellID> > faces_per_surface;
            std::map<int,std::string> surface_names;
            for(auto f_id: m_mesh.faces()){
                faces_per_surface[face_on_srf->value(f_id)].push_back(f_id);
                surface_names[face_on_srf->value(f_id)]= face_on_srf->getName();
            }
            for(auto f2s: faces_per_surface){

                FACSurface* s =
                        new FACSurface(&m_mesh,f2s.second,surface_names[f2s.first]);
                m_surfaces.push_back(s);
                m_map_face_var_to_surf[f2s.first]=s->id();


            }

            //We build a single volume!
            FACVolume* v = new FACVolume();
            m_volumes.push_back(v);
            buildTopologicalConnectivity();
        }

/*----------------------------------------------------------------------------*/
        void FACManager::initAndLinkFrom3DMesh(Mesh *AFromMesh,
                                               GeomMeshLinker *ALinker)
        {

            //the linker is updated before starting
            ALinker->clear();
            ALinker->setMesh(AFromMesh);
            ALinker->setGeometry(this);

            BoundaryExtractor3D boundary_extractor(AFromMesh, &m_mesh);

            Variable<int>* node_on_pnt = m_mesh.newVariable<int, GMDS_NODE>("on_point"  );
            Variable<int>* node_on_crv = m_mesh.newVariable<int, GMDS_NODE>("on_curve"  );
            Variable<int>* node_on_srf = m_mesh.newVariable<int, GMDS_NODE>("on_surface");
            Variable<int>* edge_on_crv = m_mesh.newVariable<int, GMDS_EDGE>("on_curve"  );
            Variable<int>* face_on_srf = m_mesh.newVariable<int, GMDS_FACE>("on_surface");

            boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, face_on_srf, node_on_srf,node_on_crv);

            //mapping from the skin mesh to the volume mesh.
            std::map<TCellID,TCellID > node_map_from_vol, edge_map_from_vol, face_map_from_vol;
            std::map<TCellID,TCellID > node_map_from_surf, edge_map_from_surf, face_map_from_surf;
            boundary_extractor.setMappings(&node_map_from_vol , &edge_map_from_vol , &face_map_from_vol,
                                           &node_map_from_surf, &edge_map_from_surf, &face_map_from_surf);
            boundary_extractor.execute();

            //to update the cad model from m_mesh, we have to add cell groups into m_mesh.

            gmds::MeshDoctor doc(&m_mesh);
            doc.buildEdgesAndX2E();
            doc.updateUpwardConnectivity();
            //DIFFERENT FROM MINE

            // vertices
            for(auto n_id: m_mesh.nodes()){
                if(node_on_pnt->value(n_id)!=0){
                    FACPoint* p = new FACPoint(&m_mesh, n_id, node_on_pnt->getName());
                    m_points.push_back(p);
                    m_map_node_var_to_point[node_on_pnt->value(n_id)]=p->id();
                }
            }

            // curves
            std::map<int,std::vector<TCellID> > edges_per_curve;
            std::map<int,std::string> curve_names;
            for(auto e_id: m_mesh.edges()) {
                if (edge_on_crv->value(e_id) != 0) {
                    edges_per_curve[edge_on_crv->value(e_id)].push_back(e_id);
                    curve_names[edge_on_crv->value(e_id)] = edge_on_crv->getName();
                }
            }
            for(auto e2c: edges_per_curve) {
                std::vector<TCellID> curve_edges = e2c.second;

                //we get all the nodes of the edge
                std::set<TCellID> all_ids;
                for (auto e_id: curve_edges){
                    Edge e = m_mesh.get<Edge>(e_id);
                    std::vector<TCellID> n_ids = e.getIDs<Node>();
                    all_ids.insert(n_ids[0]);
                    all_ids.insert(n_ids[1]);
                }
                std::vector<TCellID> node_ids;
                node_ids.assign(all_ids.begin(),all_ids.end());
                FACCurve* c =
                        new FACCurve(&m_mesh,node_ids,curve_edges, curve_names[e2c.first]);
                m_curves.push_back(c);
                m_map_edge_var_to_curve[e2c.first]=c->id();

            }

            std::map<int,std::vector<TCellID> > faces_per_surface;
            std::map<int,std::string> surface_names;
            for(auto f_id: m_mesh.faces()){
                faces_per_surface[face_on_srf->value(f_id)].push_back(f_id);
                surface_names[face_on_srf->value(f_id)]= face_on_srf->getName();
            }
            for(auto f2s: faces_per_surface){

                FACSurface* s =
                        new FACSurface(&m_mesh,f2s.second,surface_names[f2s.first]);
                m_surfaces.push_back(s);

                m_map_face_var_to_surf[f2s.first]=s->id();
                //link the volume mesh node to the new geom curve

                for(auto f_id: f2s.second) {
                    std::vector<TCellID> f_nids = m_mesh.get<Face>(f_id).getIDs<Node>();
                    for(auto n_id:f_nids) {
                        ALinker->linkNodeToSurface(node_map_from_surf[n_id], s->id());
                    }
                    std::vector<TCellID> f_eids = m_mesh.get<Face>(f_id).getIDs<Edge>();
                    for(auto e_id:f_eids) {
                        ALinker->linkEdgeToSurface(edge_map_from_surf[e_id], s->id());
                    }
                    ALinker->linkFaceToSurface(face_map_from_surf[f_id], s->id());
                }
            }

            //Now we link node on edges, then nodes on point. We do that in this way to
            //overload face linking first, then curve linking.
            for(auto c: m_curves) {
                //link the volume mesh node to the new geom curve
                std::vector<Node> cns;
                c->getMeshNodes(cns);
                for(auto n: cns) {
                    ALinker->linkNodeToCurve(node_map_from_surf[n.id()], c->id());
                }
                std::vector<Edge> ces;
                c->getMeshEdges(ces);
                for(auto e:ces) {
                    ALinker->linkEdgeToCurve(edge_map_from_surf[e.id()],c->id());
                }
            }
            for(auto p: m_points) {
                //link the volume mesh node to the new geom curve
                ALinker->linkNodeToPoint(node_map_from_surf[p->getNode().id()], p->id());

            }

            //We build a single volume!
            FACVolume* v = new FACVolume();
            m_volumes.push_back(v);
            buildTopologicalConnectivity();
        }
/*----------------------------------------------------------------------------*/
        void FACManager::buildTopologicalConnectivity(const int ADim) {

            Variable<int>* node_on_pnt = m_mesh.getVariable<int, GMDS_NODE>("on_point"  );
            Variable<int>* edge_on_crv = m_mesh.getVariable<int, GMDS_EDGE>("on_curve"  );
            Variable<int> *face_on_srf = nullptr;
            if(ADim==3) {
                face_on_srf = m_mesh.getVariable<int, GMDS_FACE>("on_surface");
            }
            FACVolume* single_vol = m_volumes[0];
            for(auto p:m_points){
                Node n = p->getNode();
                //We get the adjacent faces and edges
                std::vector<TCellID> adj_edge_ids = n.getIDs<Edge>();
                std::vector<TCellID> adj_face_ids = n.getIDs<Face>();
                // CURVE <-> POINT Connectivity
                std::set<TInt> adj_curve_ids;
                for(auto eid:adj_edge_ids){
                    if(edge_on_crv->value(eid)!=0){
                        adj_curve_ids.insert(m_map_edge_var_to_curve[edge_on_crv->value(eid)]);
                    }
                }
                for(auto cid:adj_curve_ids) {
                    GeomCurve *e = getCurve(cid);
                    e->points().push_back(p);
                    p->curves().push_back(e);
                }
                // SURFACE <-> POINT Connectivity
                std::set<TInt> adj_surface_ids;
                for(auto fid:adj_face_ids) {
                    if (ADim==3 &&face_on_srf->value(fid) != 0) {
                        adj_surface_ids.insert(m_map_face_var_to_surf[face_on_srf->value(fid)]);
                    }
                }
                for(auto sid:adj_surface_ids){
                    GeomSurface* s = getSurface(sid);
                    s->points().push_back(p);
                    p->surfaces().push_back(s);
                }
                // VOL <-> POINT Connectivity
                single_vol->points().push_back(p);
                p->volumes().push_back(single_vol);
            }

            for(auto c:m_curves){
                std::vector<Edge> c_edges;
                c->getMeshEdges(c_edges);
                // SURFACE <-> CURVE Connectivity
                std::set<TInt> adj_surface_ids;
                for(auto e:c_edges) {
                    std::vector<TCellID> adj_f_ids = e.getIDs<Face>();
                    for(auto fid:adj_f_ids) {
                        if (ADim==3 && face_on_srf->value(fid) != 0) {
                            adj_surface_ids.insert(m_map_face_var_to_surf[face_on_srf->value(fid)]);
                        }
                    }
                }
                for(auto sid:adj_surface_ids){
                    GeomSurface* s = getSurface(sid);
                    s->curves().push_back(c);
                    c->surfaces().push_back(s);
                }

                // CURVE <-> VOL Connectivity
                single_vol->curves().push_back(c);
                c->volumes().push_back(single_vol);
            }

            for(auto s:m_surfaces){
                // VOL <-> SURF Connectivity
                single_vol->surfaces().push_back(s);
                s->volumes().push_back(single_vol);
            }


        }
/*----------------------------------------------------------------------------*/
        void FACManager::initAndLinkFrom2DMesh(Mesh *AFromMesh,
                                               GeomMeshLinker *ALinker)
        {

            //the linker is updated before starting
            ALinker->clear();
            ALinker->setMesh(AFromMesh);
            ALinker->setGeometry(this);

            BoundaryExtractor2D boundary_extractor(AFromMesh, &m_mesh);
            if(!boundary_extractor.isValid()){
                throw GMDSException("FACManager::initAndLinkFrom2DMesh: "
                                    "Invalid mesh models for the 2D boundary extraction");
            }
            Variable<int>* node_on_pnt = m_mesh.newVariable<int, GMDS_NODE>("on_point"  );
            Variable<int>* node_on_crv = m_mesh.newVariable<int, GMDS_NODE>("on_curve"  );
            Variable<int>* edge_on_crv = m_mesh.newVariable<int, GMDS_EDGE>("on_curve"  );

            boundary_extractor.setColorOption(node_on_pnt, edge_on_crv, node_on_crv);

            //mapping from the skin mesh to the surface mesh.
            std::map<TCellID,TCellID > node_map_from_vol, edge_map_from_vol;
            std::map<TCellID,TCellID > node_map_from_surf, edge_map_from_surf;
            boundary_extractor.setMappings(&node_map_from_vol , &edge_map_from_vol,
                                           &node_map_from_surf, &edge_map_from_surf);
            boundary_extractor.execute();

            //to update the cad model from m_mesh, we have to add cell groups into m_mesh.

            gmds::MeshDoctor doc(&m_mesh);
            doc.buildEdgesAndX2E();
            doc.updateUpwardConnectivity();

            std::vector<TCellID> surface_nodes;
            // vertices
            for(auto n_id: m_mesh.nodes()){
                if(node_on_pnt->value(n_id)!=0){
                    FACPoint* p = new FACPoint(&m_mesh, n_id, node_on_pnt->getName());
                    m_points.push_back(p);
                    m_map_node_var_to_point[node_on_pnt->value(n_id)]=p->id();
                }
                else if(node_on_crv->value(n_id)==0){
                    //means node in surface
                    surface_nodes.push_back(n_id);
                }
            }

            // curves
            std::map<int,std::vector<TCellID> > edges_per_curve;
            std::map<int,std::string> curve_names;
            for(auto e_id: m_mesh.edges()) {
                if (edge_on_crv->value(e_id) != 0) {
                    edges_per_curve[edge_on_crv->value(e_id)].push_back(e_id);
                    curve_names[edge_on_crv->value(e_id)] = edge_on_crv->getName();
                }
            }
            for(auto e2c: edges_per_curve) {
                std::vector<TCellID> curve_edges = e2c.second;

                //we get all the nodes of the edge
                std::set<TCellID> all_ids;
                for (auto e_id: curve_edges){
                    Edge e = m_mesh.get<Edge>(e_id);
                    std::vector<TCellID> n_ids = e.getIDs<Node>();
                    all_ids.insert(n_ids[0]);
                    all_ids.insert(n_ids[1]);
                }
                std::vector<TCellID> node_ids;
                node_ids.assign(all_ids.begin(),all_ids.end());
                FACCurve* c =
                        new FACCurve(&m_mesh,node_ids,curve_edges, curve_names[e2c.first]);
                m_curves.push_back(c);
                m_map_edge_var_to_curve[e2c.first]=c->id();

            }


            //We build a single volume, even if it is a non-sense in 2D. It is
            //necessary for the buildTopo
            FACVolume* v = new FACVolume();
            m_volumes.push_back(v);
            //we also buid a single surface without more information
            FACSurface* s = new FACSurface(&m_mesh, surface_nodes, "single_surf");
            for(auto n_id: surface_nodes){
                ALinker->linkNodeToSurface(node_map_from_surf[n_id],s->id());
            }
            //Now we link node on edges, then nodes on point. We do that in this way to
            //overload face linking first, then curve linking.
            for(auto c: m_curves) {
                //link the volume mesh node to the new geom curve
                std::vector<Node> cns;
                c->getMeshNodes(cns);
                for(auto n: cns) {
                    ALinker->linkNodeToCurve(node_map_from_surf[n.id()], c->id());
                }
                std::vector<Edge> ces;
                c->getMeshEdges(ces);
                for(auto e:ces) {
                    ALinker->linkEdgeToCurve(edge_map_from_surf[e.id()],c->id());
                }
            }
            for(auto p: m_points) {
                //link the volume mesh node to the new geom curve
                ALinker->linkNodeToPoint(node_map_from_surf[p->getNode().id()], p->id());
            }

            buildTopologicalConnectivity(2);
        }

        /*----------------------------------------------------------------------------*/
        void FACManager::updateFromMesh()
        {

            std::map<TCellID,FACPoint* > map_node2point;
            std::map<TCellID,FACCurve* > map_edge2curve;

            // vertices
            for (gmds::Mesh::group_iterator<gmds::Node> itn = m_mesh.groups_begin<gmds::Node>(); itn != m_mesh.groups_end<gmds::Node>(); itn++)
            {
                gmds::CellGroup<gmds::Node>*  current_cloud = *itn;
                std::vector<TCellID> nodeIDs = current_cloud->cells();

                if(0 == nodeIDs.size()) {
                    throw GMDSException("FacetedGeomManager::importVTK a cloud is of size 0.");
                }

                // vertices are clouds with only one node
                if(1 == nodeIDs.size()) {
                    FACPoint* p = new FACPoint (&m_mesh, nodeIDs[0], current_cloud->name());
                    m_points.push_back(p);
                    map_node2point[nodeIDs[0]] = p;
                }

            }

            // in order to build our curves and surfaces, we need to rebuild some edges
            //and some connectivities
            gmds::MeshDoctor doc(&m_mesh);
            doc.buildEdgesAndX2E();
            doc.updateUpwardConnectivity();

            // WARNING numerically unstable
//	reorient();

            // curves

            // we begin by creating the surfaces in order to associate
            // faces to surfaces
            gmds::Variable<FACSurface* >* surfaceAssociation =
                    m_mesh.newVariable<FACSurface*, gmds::GMDS_FACE> ("surfaceAssociation");
            {
                gmds::Mesh::group_iterator<gmds::Face> itf = m_mesh.groups_begin<gmds::Face>();
                for (; itf != m_mesh.groups_end<gmds::Face>(); itf++)
                {
                    gmds::CellGroup<gmds::Face>*  current_surface = *itf;
                    std::vector<TCellID> surf_faces = current_surface->cells();

                    if(0 == surf_faces.size()) {
                        throw GMDSException("FacetedGeomManager::importVTK a surface is of size 0.");
                    }

                    std::vector<FACPoint* > points;
                    std::vector<FACCurve* > curves;
                    FACSurface* s =
                            new FACSurface(&m_mesh, surf_faces, current_surface->name());
                    m_surfaces.push_back(s);

                    for(unsigned int iFace=0; iFace<surf_faces.size(); iFace++) {
                        (*surfaceAssociation)[surf_faces[iFace]] = s;
                    }
                } // for(; its != mesh_.surfaces_end(); its++) {

            }

            // now we get create the curves
            for (gmds::Mesh::group_iterator<gmds::Node> itn = m_mesh.groups_begin<gmds::Node>(); itn != m_mesh.groups_end<gmds::Node>(); itn++)
            {
                gmds::CellGroup<gmds::Node>*  current_cloud = *itn;
                std::vector<TCellID> nodeIDs = current_cloud->cells();

                // curves are clouds with more than one node
                if(1 != nodeIDs.size()) {

                    std::vector<TCellID> orderedNodeIDs;

                    // we need to order the nodes
                    // first we find the two end-nodes
                    TCellID firstNodeID = NullID;
                    TCellID lastNodeID  = NullID;

                    for(unsigned int iNode=0; iNode<nodeIDs.size(); iNode++) {
                        if(map_node2point.end() != map_node2point.find(nodeIDs[iNode])) {
                            if(NullID == firstNodeID) {
                                firstNodeID = nodeIDs[iNode];
                            }
                            else
                            {
                                lastNodeID = nodeIDs[iNode];
                            }
                        }
                    }
//			std::cout<<"Curve from "<<firstNodeID<<" to "<<lastNodeID<<std::endl;

                    if(NullID == firstNodeID) {
                        throw GMDSException(""
                                            "FacetedGeomManager::importVTK point of firstNode of curve was not found.");
                    }
                    // in this case we have a loop curve
                    if(NullID == lastNodeID) {
                        lastNodeID = firstNodeID;
                    }

                    // mark all the nodes of this curve
                    gmds::Variable<int>* markCurveNodes = m_mesh.newVariable<int, gmds::GMDS_NODE>("MarkCurveNodes");

                    for(unsigned int iNode=0; iNode<nodeIDs.size(); iNode++) {
                        markCurveNodes->set(nodeIDs[iNode],1);
                    }

                    // identify if this curve is adjacent to only one surface;
                    bool hasOneSurface = true;
                    // for example in a cylinder
                    for(unsigned int iNode=0; iNode<nodeIDs.size(); iNode++) {

                        std::vector<Edge> edges = m_mesh.get<gmds::Node>(nodeIDs[iNode]).get<Edge>();

                        for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
                            std::vector<Node> edgeNodes = edges[iEdge].get<Node>();
                            if((*markCurveNodes)[edgeNodes[0].id()]==1 &&
                               (*markCurveNodes)[edgeNodes[1].id()]==1)
                            {
                                std::vector<Face> edgeFaces = edges[iEdge].get<Face>();
                                if(2 < edgeFaces.size()) {
                                    throw GMDSException("FacetedGeomManager::importVTK "
                                                        "edge must have two faces");
                                }
                                else if(1 == edgeFaces.size()){
                                    hasOneSurface=true;
                                }
                                else if(2==edgeFaces.size()) {
                                    if ((*surfaceAssociation)[edgeFaces[0].id()] !=
                                        (*surfaceAssociation)[edgeFaces[1].id()]) {
                                        hasOneSurface = false;
                                        break;
                                    }
                                }
                            }
                        }
                        if(!hasOneSurface) {
                            break;
                        }
                    }

                    std::vector<gmds::TCellID> edgeIDs_path;
                    TCellID current_edge_id = NullID;

                    std::vector<gmds::TCellID> orderedNodes;
                    TCellID current_node_id = firstNodeID;
                    gmds::Node  current_node = m_mesh.get<Node>(current_node_id);
                    orderedNodes.push_back(current_node.id());

                    // we build the path of edges
                    // we use a do...while here because firstNode can be equal
                    // to lastNode in case of a loop curve
                    do {
                        std::vector<Edge> edges = current_node.get<Edge>();

                        // select next edge
                        bool edgeFound = false;
                        for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {

                            // must be different from the current one
                            if(edges[iEdge].id() != current_edge_id) {

                                // its two nodes must be on the curve
                                std::vector<Node> edgeNodes = edges[iEdge].get<Node>();
                                if((*markCurveNodes)[edgeNodes[0].id()]==1 &&
                                   (*markCurveNodes)[edgeNodes[1].id()]==1)
                                {
                                    // this must be a curve edge
                                    std::vector<Face> edgeFaces = edges[iEdge].get<Face>();

                                    if(2 < edgeFaces.size()) {
                                        throw GMDSException("FacetedGeomManager::importVTK "
                                                            "edge must have two faces");
                                    }

                                    // in case of a one-surface edge, no need for this test
                                    if(	hasOneSurface ||
                                           ((*surfaceAssociation)[edgeFaces[0].id()] != (*surfaceAssociation)[edgeFaces[1].id()])) {

                                        edgeFound = true;
                                        Edge current_edge = edges[iEdge];
                                        edgeIDs_path.push_back(current_edge.id());
                                        current_edge_id = current_edge.id();

                                        if(edgeNodes[0] == current_node) {
                                            current_node = edgeNodes[1];
                                        } else {
                                            current_node = edgeNodes[0];
                                        }
                                        current_node_id = current_node.id();
                                        orderedNodes.push_back(current_node.id());
                                        break;
                                    }
                                }
                            }
                        } // for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {

                        if(!edgeFound) {
                            throw GMDSException("FacetedGeomManager::importVTK "
                                                "next edge not found.");
                        }

                    } while (current_node_id != lastNodeID);

                    m_mesh.deleteVariable(GMDS_NODE,"MarkCurveNodes");

                    FACPoint* p1 = map_node2point[firstNodeID];
                    FACPoint* p2 = map_node2point[lastNodeID];

                    FACCurve* c =
                            new FACCurve(&m_mesh, p1, p2, orderedNodes, edgeIDs_path, current_cloud->name());

                    m_curves.push_back(c);
//                    p1->add(c);
//                    p2->add(c);

                    for(unsigned int iEdge=0; iEdge<edgeIDs_path.size(); iEdge++) {
                        map_edge2curve[edgeIDs_path[iEdge]] = c;
                    }

                } // if(1 != nodes.size()) {

            } // curves

//            // surfaces
//            // they were previously created;
//            // now we will add adjacency relations
//            for(unsigned int iSurf=0; iSurf<m_surfaces.size(); iSurf++) {
//                std::vector<Face> surf_faces;
//                m_surfaces[iSurf]->getMeshFaces(surf_faces);
//
//                if(0 == surf_faces.size()) {
//                    throw GMDSException("FacetedGeomManager::importVTK a surface is of size 0.");
//                }
//
//                // we get all the curves surrounding this surfaces
//                // we get the vertices as well
//                std::set<FacetedCurve*> curves_set;
//                std::vector<FacetedCurve*> curves;
//                std::set<FacetedPoint*> points_set;
//                std::vector<FacetedPoint*> points;
//
//                for(unsigned int iFace=0; iFace<surf_faces.size(); iFace++) {
//                    std::vector<Edge> edges = surf_faces[iFace].get<Edge>();
//
//                    for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
//                        if(map_edge2curve.end() != map_edge2curve.find(edges[iEdge].getID())) {
//
//                            FacetedCurve* curv = map_edge2curve[edges[iEdge].getID()];
//                            curves_set.insert(curv);
//                            if(NULL != curv->getFirstPoint()) {
//                                points_set.insert(dynamic_cast<gmds::geom::FacetedPoint* > (curv->getFirstPoint()));
//                            }
//                            if(NULL != curv->getSecondPoint()) {
//                                points_set.insert(dynamic_cast<gmds::geom::FacetedPoint* > (curv->getSecondPoint()));
//                            }
//                        }
//                    }
//                }
//
//                curves.insert(curves.begin(),curves_set.begin(),curves_set.end());
//                points.insert(points.begin(),points_set.begin(),points_set.end());
//
//                m_surfaces[iSurf]->replace(points);
//                m_surfaces[iSurf]->replace(curves);
//
//                //connectivity curves -> surfaces
//                for(unsigned int i=0; i<curves.size(); i++) {
//                    curves[i]->add(m_surfaces[iSurf]);
//                }
//
//            } // for(unsigned int iSurf=0; iSurf<m_surfaces.size(); iSurf++) {

            m_mesh.deleteVariable(GMDS_FACE,"surfaceAssociation");

//            FacetedVolume* v = new FacetedVolume(m_surfaces);
//            m_volumes.push_back(v);

//            // connectivity surface -> volume
//            for(unsigned int i=0; i<m_surfaces.size(); i++) {
//                m_surfaces[i]->add(v);
//            }

            // we delete every clouds and surfaces of mesh_

//            // vertices
//            {
//                unsigned int nbClouds = mesh_.getNbClouds();
//                for(unsigned int iCloud=0; iCloud<nbClouds; iCloud++) {
//
//                    IGMesh::clouds_iterator itc = mesh_.clouds_begin();
//                    IGMesh::cloud cl = *itc;
//                    mesh_.deleteCloud(cl);
//                }
//            }
//
//            // surfaces
//            {
//                unsigned int nbSurfaces = mesh_.getNbSurfaces();
//                for(unsigned int iSurf=0; iSurf<nbSurfaces; iSurf++) {
//
//                    IGMesh::surfaces_iterator its = mesh_.surfaces_begin();
//                    IGMesh::surface surf = *its;
//                    mesh_.deleteSurface(surf);
//                }
//            }

#ifdef _DEBUG_	 //edge [0, 1] gives 3

            std::cout<<"We imported "<<std::endl;
	        std::cout<<"nb vertices "<<this->getNbPoints()<<std::endl;
	        std::cout<<"nb curves "<<this->getNbCurves()<<std::endl;
	        std::cout<<"nb surfaces "<<this->getNbSurfaces()<<std::endl;
#endif //_DEBUG_

        }
        /*----------------------------------------------------------------------------*/
        GeomVolume* FACManager::newVolume()
        {
            FACVolume* v = new FACVolume();
            m_volumes.push_back(v);
            return v;
        }
/*----------------------------------------------------------------------------*/
        GeomSurface* FACManager::newSurface()
        {
            FACSurface* s = new FACSurface(&m_mesh);
            m_surfaces.push_back(s);
            return s;
        }
/*----------------------------------------------------------------------------*/
        GeomCurve* FACManager::newCurve()
        {
            FACCurve* c = new FACCurve(&m_mesh);
            m_curves.push_back(c);
            return c;
        }
/*----------------------------------------------------------------------------*/
        GeomPoint* FACManager::newPoint()
        {
            FACPoint* p = new FACPoint(&m_mesh);
            m_points.push_back(p);
            return p;
        }
/*----------------------------------------------------------------------------*/
        Mesh& FACManager::getMeshView()
        {
            return m_mesh;
        }
/*----------------------------------------------------------------------------*/
        TInt FACManager::getNbPoints() const
        {
            return m_points.size();
        }
/*----------------------------------------------------------------------------*/
        TInt FACManager::getNbCurves() const
        {
            return m_curves.size();
        }
/*----------------------------------------------------------------------------*/
        TInt FACManager::getNbSurfaces() const
        {
            return m_surfaces.size();
        }
/*----------------------------------------------------------------------------*/
        TInt FACManager::getNbVolumes() const
        {
            return m_volumes.size();
        }
/*----------------------------------------------------------------------------*/
        void FACManager::getVolumes(std::vector<GeomVolume*>& volumes) const
        {
            volumes.clear();
            volumes.resize(m_volumes.size());
            for(unsigned int i=0; i < m_volumes.size(); i++)
                volumes[i]=m_volumes[i];
        }
/*----------------------------------------------------------------------------*/
        void FACManager::getSurfaces(std::vector<GeomSurface*>& surfaces) const
        {
            surfaces.clear();
            surfaces.resize(m_surfaces.size());
            for(unsigned int i=0; i < m_surfaces.size(); i++)
                surfaces[i]=m_surfaces[i];
        }
/*----------------------------------------------------------------------------*/
        void FACManager::getCurves(std::vector<GeomCurve*>& curves) const
        {
            curves.clear();
            curves.resize(m_curves.size());
            for(unsigned int i=0; i < m_curves.size(); i++)
                curves[i]=m_curves[i];
        }
/*----------------------------------------------------------------------------*/
        void FACManager::
        getPoints(std::vector<GeomPoint*>& points) const
        {
            points.clear();
            points.resize(m_points.size());
            for(unsigned int i=0; i < m_points.size(); i++)
                points[i]=m_points[i];
        }
/*----------------------------------------------------------------------------*/
        GeomPoint* FACManager::getPoint(const gmds::TInt AID) {
            for(auto p:m_points){
                if(p->id()==AID)
                    return p;
            }
            return nullptr;
        }
/*----------------------------------------------------------------------------*/
        GeomCurve* FACManager::getCurve(const gmds::TInt AID) {
            for(auto c:m_curves){
                if(c->id()==AID)
                    return c;
            }
            return nullptr;
        }
/*----------------------------------------------------------------------------*/
        GeomSurface* FACManager::getSurface(const gmds::TInt AID) {
            for(auto s:m_surfaces){
                if(s->id()==AID)
                    return s;
            }
            return nullptr;}
/*----------------------------------------------------------------------------*/
        GeomVolume* FACManager::getVolume(const gmds::TInt AID) {
            for(auto v:m_volumes){
                if(v->id()==AID)
                    return v;
            }
            return nullptr;}
/*----------------------------------------------------------------------------*/
        int FACManager::getCommonCurve(GeomPoint *AP1,
                                       GeomPoint *AP2) const {

            FACPoint* p1 = dynamic_cast<FACPoint*>(AP1);
            FACPoint* p2 = dynamic_cast<FACPoint*>(AP2);

            Node n1 = p1->getNode();
            Node n2 = p2->getNode();

            Variable<int>* edge_on_crv = m_mesh.getVariable<int, GMDS_EDGE>("on_curve"  );

            std::vector<Edge> e1 = n1.get<Edge>();
            std::vector<Edge> e2 = n2.get<Edge>();

            for(auto e1_edge:e1){
                for(auto e2_edge:e2){
                    if((*edge_on_crv)[e1_edge.id()]==(*edge_on_crv)[e2_edge.id()]){
                        return (*edge_on_crv)[e1_edge.id()];
                    }
                }
            }
            return -1;
        }
/*----------------------------------------------------------------------------*/
        int FACManager::getCommonSurface(GeomCurve *AC1,
                                         GeomCurve *AC2) const {
            FACCurve* c1 = dynamic_cast<FACCurve*>(AC1);
            FACCurve* c2 = dynamic_cast<FACCurve*>(AC2);
            std::vector<Node> n1, n2;
            c1->getMeshNodes(n1);
            c2->getMeshNodes(n2);

            Variable<int>* node_on_srf = m_mesh.getVariable<int, GMDS_NODE>("on_surface");

            std::set<int> surf_1;
            for(auto n1_node:n1){
                std::vector<Edge> n1_edge = n1_node.get<Edge>();
                for(auto e:n1_edge){
                    std::vector<TCellID > e_nodes = e.getIDs<Node>();
                    TCellID opp_node=(e_nodes[0]==n1_node.id())?e_nodes[1]:e_nodes[0];
                    int surf_color = (*node_on_srf)[opp_node];
                    if(surf_color!=0){
                        surf_1.insert(surf_color);
                    }
                }
            }
            std::set<int> surf_2;
            for(auto n2_node:n2){
                std::vector<Edge> n2_edge = n2_node.get<Edge>();
                for(auto e:n2_edge){
                    std::vector<TCellID > e_nodes = e.getIDs<Node>();
                    TCellID opp_node=(e_nodes[0]==n2_node.id())?e_nodes[1]:e_nodes[0];
                    int surf_color = (*node_on_srf)[opp_node];
                    if(surf_color!=0){
                        surf_2.insert(surf_color);
                    }
                }
            }

            //surf_1 contains all the surf id around curve 1 (idem for surf_2)
            for(auto s1:surf_1){
                for(auto s2:surf_2){
                    if(s1==s2)
                        return s1;
                }
            }
            return -1;
        }
        /*----------------------------------------------------------------------------*/
        void
        FACManager::write_surfaces(std::string AFilename) const
        {
            // build the temporary mesh
            gmds::Mesh mesh_tmp(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
            std::map<gmds::TCellID, gmds::TCellID> old2NewNodes;
            std::map<gmds::TCellID, gmds::TCellID> old2NewFaces;

            for(auto id: m_mesh.nodes()){
                gmds::Node n = m_mesh.get<gmds::Node>(id);

                gmds::Node  n_new = mesh_tmp.newNode(n.X(), n.Y(), n.Z());
                old2NewNodes.emplace(n.id(), n_new.id());
            }

            for(auto id: m_mesh.faces()) {
                gmds::Face f = m_mesh.get<gmds::Face>(id);

                std::vector<gmds::TCellID> nodeIDs = f.getIDs<gmds::Node>();

                for(int i_n=0; i_n<nodeIDs.size(); i_n++) {
                    nodeIDs[i_n] = old2NewNodes[nodeIDs[i_n]];
                }

                gmds::Face f_new = mesh_tmp.newFace(nodeIDs);
                old2NewFaces.emplace(f.id(), f_new.id());
            }

            // create variable
            gmds::Variable<int> *v = mesh_tmp.newVariable<int, GMDS_FACE>("material_interface");
            gmds::Variable<math::Vector3d> *n = mesh_tmp.newVariable<math::Vector3d, GMDS_FACE>("normal");

            // and fill it
            for(unsigned int iSurf=0; iSurf < m_surfaces.size(); iSurf++) {
                std::vector<Face> surf_faces;
                m_surfaces[iSurf]->getMeshFaces(surf_faces);

                for (auto f: surf_faces) {
                    (*v)[old2NewFaces[f.id()]] = iSurf;
                    n->set(f.id(),f.normal());
                }
            }

            gmds::IGMeshIOService ioService(&mesh_tmp);
            gmds::VTKWriter vtkWriter(&ioService);
            vtkWriter.setCellOptions(gmds::N|gmds::F);
            vtkWriter.setDataOptions(gmds::N|gmds::F);
            vtkWriter.write(AFilename);
        }
/*----------------------------------------------------------------------------*/
        FACSurface* FACManager::getFACSurface(const gmds::TInt AID){
            for(auto s:m_surfaces){
                if(s->id()==AID)
                    return s;
            }
            return nullptr;
        }
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
