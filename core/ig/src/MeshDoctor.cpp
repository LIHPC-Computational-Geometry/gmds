/*----------------------------------------------------------------------------*/
/*
 * MeshDoctor.cpp
 *
 *  Created on: 22 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/MeshDoctor.h>
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    MeshDoctor::MeshDoctor(Mesh* AMesh)
    : m_mesh(AMesh)
    {}
    /*----------------------------------------------------------------------------*/
    MeshDoctor::~MeshDoctor()
    = default;
    
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::setMesh(Mesh* AMesh)
    {
        m_mesh=AMesh;
    }
    /*----------------------------------------------------------------------------*/
    TCoord MeshDoctor::isLeft(Node& AN1, Node& AN2, Node& AN3)
    {
        return ( (AN2.X()-AN1.X()) * (AN3.Y()-AN1.Y()) -
                (AN3.X()-AN1.X()) * (AN2.Y()-AN1.Y()) );
    }
    /*----------------------------------------------------------------------------*/
    int MeshDoctor::orient2DFaces()
    {
        
        TInt nb_reorientation =0;
        for(auto f_id:m_mesh->faces())
        {
            Face f = m_mesh->get<Face>(f_id);
            if (orient2DFace(f))
                nb_reorientation++;
        }
        
        return nb_reorientation;
    }
    /*------------------------------------------------------------------------*/
    bool MeshDoctor::orient2DFace(Face& AF)
    {
	     bool isReoriented = false;
	     std::vector<Node> nodes = AF.get<Node>();
	     TCoord orientation=0;
	     if(AF.type()==GMDS_TRIANGLE) {
		     orientation = isLeft(nodes[0],nodes[1],nodes[2]);
	     }
	    else {
		     //find the rightmost lowest vertex of the polygon
		     unsigned int index_min=0;
		     TCoord x_min = nodes[0].X();
		     TCoord y_min = nodes[0].Y();
		     for(unsigned int
		             i=0;i<nodes.size();i++)
		     {
			     if(nodes[i].Y()>y_min)
				     continue;
			     if(nodes[i].Y()==y_min) {	// just as low
				     if(nodes[i].X()<x_min)  // and to left
					     continue;
			     }

			     index_min =i;
			     x_min = nodes[i].X();
			     y_min = nodes[i].Y();
		     }

		     if(index_min==0)
			     orientation = isLeft(nodes[nodes.size()-1],nodes[0],nodes[1]);
		     else if (index_min==nodes.size()-1)
			     orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[0]);
		     else
			     orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[index_min+1]);
	     }
	     if(orientation>0.0) // clockwise or degenerated (=0)
	     {
		     isReoriented= true;
		     std::vector<Node> nodes_inv;
		     nodes_inv.resize(nodes.size());
		     auto node_size = nodes.size();
		     for(unsigned int i=0;i<node_size;i++)
			     nodes_inv[i] = nodes[node_size-1-i];
		     AF.set<Node>(nodes_inv);
	     }

	     return isReoriented;
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildFacesAndR2F() const
    {
        if(m_mesh->getModel().has(F) && m_mesh->getModel().has(R2F)) {
            buildFAndR2F();
        } else {
            throw GMDSException("MeshDoctor::buildFacesAndR2F "
                                "Can not be called when mesh model does not have F && R2F");
        }
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildEdgesAndX2E() const
    {
        if(m_mesh->getModel().has(E) && (m_mesh->getModel().has(R2E) || m_mesh->getModel().has(F2E))) {
            buildE();
            if(m_mesh->getModel().has(R2E)) {
                buildR2E(m_mesh->getModel());
            }
            if(m_mesh->getModel().has(F2E)) {
                buildF2E(m_mesh->getModel());
            }
        } else {
            throw GMDSException("MeshDoctor::buildEdgesAndX2E "
                                "Can not be called when mesh model does not have E && (R2E || F2E)");
        }
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildBoundaryCells() const
    {
        //TODO: define how to implement it
        // To build boundary edge, we need to find boundary nodes in each cell
        // we require so the N2F or N2R connectivity
        if(m_mesh->getModel().has(E) &&
           m_mesh->getModel().has(F) &&
           m_mesh->getModel().has(R))  {
            //3D CASE
            throw GMDSException("MeshDoctor::buildBoundaryCells "
                                "3D Case not yet implemented");
        } else if(m_mesh->getModel().has(E) &&
                  m_mesh->getModel().has(F) )  {
            //2D CASE
            if(!m_mesh->getModel().has(N2F)){
                throw GMDSException("MeshDoctor::buildBoundaryCells "
                                    "2D case require N2F connection");
            }
            //We build boundary edges only
            for(auto f_id:m_mesh->faces()){
                Face f = m_mesh->get<Face>(f_id);
                std::vector<Node> ns = f.get<Node>();
                for(auto i=0; i<ns.size();i++){
                    Node ni=ns[i];
                    Node nj=ns[(i+1)%ns.size()];
                    std::vector<TCellID> fij = m_mesh->getCommonFaces(ni,nj);
                    if(fij.size()==1){
                        //Means we are on a boundary edge, we build it
                        m_mesh->newEdge(ni,nj);
                    }
                }
            }
        } else  {
            throw GMDSException("MeshDoctor::buildBoundaryCells "
                                "Can not be used when mesh model does not have E, F, and R in 3D or E and F in 2D");
        }
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::updateUpwardConnectivity() const
    {
        MeshModel mod = m_mesh->getModel();
        if(mod.has(N2F) && mod.has(F2N))
        {
            for(auto id_f : m_mesh->faces()) {
                Face f = m_mesh->get<Face>(id_f);
                std::vector<Node> nodes = f.get<Node>();
                for(auto & node : nodes)
                    if(!node.has<Face>(f)) {
                        node.add<Face>(f);
                    }
            }
            
        }
        if(mod.has(E2F) && mod.has(F2E))
        {
            for(auto id_f : m_mesh->faces()) {
                Face f = m_mesh->get<Face>(id_f);
                std::vector<Edge> edges = f.get<Edge>();
                
                for(auto & edge : edges)
                    if(!edge.has<Face>(f)) {
                        edge.add<Face>(f);
                    }
            }
        }
        if(mod.has(N2R) && mod.has(R2N))
        {

            for(auto id_r : m_mesh->regions()) {
                Region r = m_mesh->get<Region>(id_r);
                std::vector<Node> nodes = r.get<Node>();
                
                for(auto & node : nodes)
                    if(!node.has<Region>(r)) {
                        node.add<Region>(r);
                    }
            }
        }
        
        if(mod.has(E2R) && mod.has(R2E))
        {
            for(auto id_r : m_mesh->regions()) {
                Region r = m_mesh->get<Region>(id_r);
                std::vector<Edge> edges = r.get<Edge>();
                
                for(auto & edge : edges)
                    if(!edge.has<Region>(r)) {
                        edge.add<Region>(r);
                    }
            }
        }
        if(mod.has(N2E) && mod.has(E2N))
        {
            for(auto e_id : m_mesh->edges())
            {
                Edge e = m_mesh->get<Edge>(e_id);
                std::vector<Node> nodes = e.get<Node>();
                
                for(auto & node : nodes)
                    if(!node.has<Edge>(e)) {
                        node.add<Edge>(e);
                    }
            }
        }
        
        if(mod.has(F2R) && mod.has(R2F))
        {
            for(auto id_r : m_mesh->regions()) {
                Region r = m_mesh->get<Region>(id_r);
                std::vector<Face> faces = r.get<Face>();
                
                for(auto & face : faces) {
                    if(!face.has<Region>(r)) {
                        face.add<Region>(r);
                    }
                }
            }
        }
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildE() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<VirtualEdge::EdgeID, TCellID> tmp_edges;
        
        // check model validity
        if(!(mod.has(E) && mod.has(E2N))) {
            throw GMDSException("MeshDoctor::buildE E or E2N missing.");
        }
        
        if(mod.has(R) && mod.has(R2N)){
            buildEfromR();
        }
        else if(mod.has(F) && mod.has(F2N)){
		      buildEfromF();
        }
        else
            throw GMDSException("MeshDoctor::buildE R or F missing");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildEfromF() const
    {
	     MeshModel mod = m_mesh->getModel();

	     std::map<VirtualEdge::EdgeID, TCellID> tmp_edges;

	     // first register the existing edges
	     if(mod.has(E) && mod.has(E2N)) {
		      for(auto id_e : m_mesh->edges())
		      {
			       Edge e = m_mesh->get<Edge>(id_e);
			       std::vector<TCellID> e_nodes = e.getIDs<Node>();
			       VirtualEdge fke(e_nodes[0], e_nodes[1]);
			       tmp_edges[fke.getID()] = e.id();
		      }
	     } else {
		      throw GMDSException("MeshDoctor::buildEfromF E or E2N missing.");
	     }

	     if(mod.has(F) && mod.has(F2N)){
		      for(auto f_id: m_mesh->faces()){
			       Face f = m_mesh->get<Face>(f_id);
			       ECellType f_type = f.type();
			       std::vector<TCellID> f_nodes = f.getIDs<Node>();
			       if(f_type==GMDS_TRIANGLE){
				        addEdge(f_nodes[0],f_nodes[1],tmp_edges);
				        addEdge(f_nodes[1],f_nodes[2],tmp_edges);
				        addEdge(f_nodes[2],f_nodes[0],tmp_edges);
			       }
			       else if(f_type==GMDS_QUAD){
				        addEdge(f_nodes[0],f_nodes[1],tmp_edges);
				        addEdge(f_nodes[1],f_nodes[2],tmp_edges);
				        addEdge(f_nodes[2],f_nodes[3],tmp_edges);
				        addEdge(f_nodes[3],f_nodes[0],tmp_edges);
			       }
			       else if(f_type==GMDS_POLYGON){
				        unsigned int nb_nodes = f_nodes.size();
				        for(unsigned int i_node=0;i_node<nb_nodes;i_node++){
					         addEdge(f_nodes[i_node],f_nodes[(i_node+1)%nb_nodes],tmp_edges);
				        }
			       }
			       else
				        throw GMDSException("MeshDoctor::buildEfromF Cell type unknown.");
		      }
	     }
	     else
		      throw GMDSException("MeshDoctor::buildEfromF F or F2N missing");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildEfromR() const
    {
	     MeshModel mod = m_mesh->getModel();

	     std::map<VirtualEdge::EdgeID, TCellID> tmp_edges;

	     // first register the existing edges
	     if(mod.has(E) && mod.has(E2N)) {
		      for(auto id_e : m_mesh->edges())
		      {
			       Edge e = m_mesh->get<Edge>(id_e);
			       std::vector<TCellID> e_nodes = e.getIDs<Node>();
			       VirtualEdge fke(e_nodes[0], e_nodes[1]);
			       tmp_edges[fke.getID()] = e.id();
		      }
	     } else {
		      throw GMDSException("MeshDoctor::buildEfromR E or E2N missing.");
	     }

	     if(mod.has(R) && mod.has(R2N)){
		      for(auto r_id : m_mesh->regions()){
			       Region r = m_mesh->get<Region>(r_id);
			       ECellType r_type = r.type();
			       std::vector<TCellID> r_nodes = r.getIDs<Node>();
			       if(r_type==GMDS_TETRA){
				        addEdge(r_nodes[0],r_nodes[1],tmp_edges);
				        addEdge(r_nodes[0],r_nodes[2],tmp_edges);
				        addEdge(r_nodes[0],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[2],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[3],tmp_edges);
			       }
			       else if(r_type==GMDS_HEX){
				        addEdge(r_nodes[0],r_nodes[1],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[2],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[3],r_nodes[0],tmp_edges);
				        addEdge(r_nodes[4],r_nodes[5],tmp_edges);
				        addEdge(r_nodes[5],r_nodes[6],tmp_edges);
				        addEdge(r_nodes[6],r_nodes[7],tmp_edges);
				        addEdge(r_nodes[7],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[0],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[5],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[6],tmp_edges);
				        addEdge(r_nodes[3],r_nodes[7],tmp_edges);
			       }
			       else if(r_type==GMDS_PYRAMID){
				        addEdge(r_nodes[0],r_nodes[1],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[2],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[3],r_nodes[0],tmp_edges);
				        addEdge(r_nodes[0],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[3],r_nodes[4],tmp_edges);
			       }
			       else if(r_type==GMDS_PRISM3){
				        addEdge(r_nodes[0],r_nodes[1],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[2],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[0],tmp_edges);
				        addEdge(r_nodes[3],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[4],r_nodes[5],tmp_edges);
				        addEdge(r_nodes[5],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[0],r_nodes[3],tmp_edges);
				        addEdge(r_nodes[1],r_nodes[4],tmp_edges);
				        addEdge(r_nodes[2],r_nodes[5],tmp_edges);
			       }
			       else
				        throw GMDSException("MeshDoctor::buildEfromR Cell type unknown.");
		      }
	     }
	     else
		      throw GMDSException("MeshDoctor::buildEfromR R or R2N missing");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::addEdge(TCellID AN1, TCellID AN2,
                               std::map<VirtualEdge::EdgeID, TCellID>& AFakeEdgeMap) const
    {
        VirtualEdge fke(AN1, AN2);
        if(AFakeEdgeMap.find(fke.getID())==AFakeEdgeMap.end())
        {
            //New face to add
            Edge e = m_mesh->newEdge(AN1,AN2);
            AFakeEdgeMap[fke.getID()]=e.id();
        }
    }
    /*----------------------------------------------------------------------------*/
    TCellID MeshDoctor::
    addFace(std::vector<TCellID>& ANodeIDs,
            std::map<VirtualFace::FaceID, TCellID>& AFakeFaceMap) const
    {
        TCellID f_id;

        VirtualFace fkf(ANodeIDs);
        auto itf =AFakeFaceMap.find(fkf.getID());

        if(itf==AFakeFaceMap.end()){
            //New face to add
            Face f = m_mesh->newFace(ANodeIDs);
            AFakeFaceMap[fkf.getID()]=f.id();
            f_id=f.id();
        }
        else{
            f_id= itf->second;
        }
        return f_id;
    }
    /*----------------------------------------------------------------------------*/
    TCellID MeshDoctor::
    addFace(Face& AFace,
            std::map<VirtualFace::FaceID, TCellID>& AFakeFaceMap)
    {
        TCellID f_id;
        std::vector<TCellID> f_nodes = AFace.getIDs<Node>();
        
        VirtualFace fkf(f_nodes);
        auto itf =AFakeFaceMap.find(fkf.getID());
        if(itf==AFakeFaceMap.end())
        {
            //New face to add
            AFakeFaceMap[fkf.getID()]=AFace.id();
            f_id = AFace.id();
        }
        else{
            f_id = itf->second;
        }
        return f_id;
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildFAndR2F() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<VirtualFace::FaceID, TCellID> tmp_faces;
        //==============================================================
        // First we put existing faces into our tmp_faces container
        //==============================================================
        if(mod.has(F) && mod.has(F2N)){

            for(auto f_id : m_mesh->faces()){
                Face f = m_mesh->get<Face>(f_id);
                std::vector<TCellID> f_nodes = f.getIDs<Node>();
                addFace(f,tmp_faces);
            }
        }
        //==============================================================
        // Second we look for missing faces
        //==============================================================
        if(mod.has(R) && mod.has(R2N)){

            for(auto r_id : m_mesh->regions()){
                Region r = m_mesh->get<Region>(r_id);
                ECellType r_type = r.type();
                std::vector<TCellID> r_nodes = r.getIDs<Node>();
                std::vector<TCellID> r_faces;
                if(r_type==GMDS_TETRA){
                    r_faces.resize(4);
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(3);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[3];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[3];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_HEX){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(6);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[5];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[6];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[3];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                    //FACE 6
                    f_nodes[0] = r_nodes[4];
                    f_nodes[1] = r_nodes[5];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[7];
                    r_faces[5] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PRISM3){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(5);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    f_nodes[3] = r_nodes[3];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    f_nodes[3] = r_nodes[5];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[5];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PYRAMID){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(5);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[4];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[4];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[4];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                }
                else
                    throw GMDSException("MeshDoctor::buildF Not yet implemented");
                r.set<Face>(r_faces);
            }//for(;!it.isDone();it.next())
        }//if(mod.has(R) && mod.has(R2N))
        else
            throw GMDSException("MeshDoctor::buildF Not yet implemented");
    }    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildF() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<VirtualFace::FaceID, TCellID> tmp_faces;
        //==============================================================
        // First we put existing faces into our tmp_faces container
        //==============================================================
        if(mod.has(F) && mod.has(F2N)){

            for(auto f_id : m_mesh->faces()){
                Face f = m_mesh->get<Face>(f_id);
                std::vector<TCellID> f_nodes = f.getIDs<Node>();
                addFace(f,tmp_faces);
            }
        }
        //==============================================================
        // Second we look for missing faces
        //==============================================================
        if(mod.has(R) && mod.has(R2N)){
            for(auto r_id : m_mesh->regions()){
                Region r = m_mesh->get<Region>(r_id);
                ECellType r_type = r.type();
                std::vector<TCellID> r_nodes = r.getIDs<Node>();
                if(r_type==GMDS_TETRA){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(3);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_HEX){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[6];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    //FACE 6
                    f_nodes[0] = r_nodes[4];
                    f_nodes[1] = r_nodes[5];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[7];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PRISM3){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    f_nodes[3] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    f_nodes[3] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PYRAMID){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                }
                else
                    throw GMDSException("MeshDoctor::buildF Not yet implemented");
            }//for(;!it.isDone();it.next())
        }//if(mod.has(R) && mod.has(R2N))
        else
            throw GMDSException("MeshDoctor::buildF Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildN2E(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(E2N)){

            for(auto e_id : m_mesh->edges()){
                Edge e = m_mesh->get<Edge>(e_id);
                std::vector<Node> e_nodes;
                e.get<Node>(e_nodes);
                for(auto n_i : e_nodes) {
                    n_i.add<Edge>(e);
                }
            }
        }
        else
            throw GMDSException("MeshDoctor::buildN2E Not yet implemented");}
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildN2F(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(F2N)){

            for(auto f_id : m_mesh->faces()){
                Face f = m_mesh->get<Face>(f_id);
                std::vector<Node> f_nodes;
                f.get<Node>(f_nodes);
                for(auto n_i : f_nodes) {
                    n_i.add<Face>(f);
                }
            }
        }
        else
            throw GMDSException("MeshDoctor::buildN2F Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildN2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2N) && ARefModel.has(N2R)){
            for(auto r_id : m_mesh->regions()){
                Region r = m_mesh->get<Region>(r_id);
                std::vector<Node> r_nodes;
                r.get<Node>(r_nodes);
                for(auto n_i : r_nodes){
                    n_i.add<Region>(r);
                }
            }
        }
        else
            throw GMDSException("MeshDoctor::buildN2R Not yet implemented for this model");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildF2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2F) && ARefModel.has(F2R)){
            for(auto r_id : m_mesh->regions()){
                Region r = m_mesh->get<Region>(r_id);
                std::vector<Face> r_faces;
                r.get<Face>(r_faces);
                for(auto f_i:r_faces){
                    f_i.add<Region>(r);
                }
            }
        }
        else
            throw GMDSException("MeshDoctor::buildF2R Not yet implemented for this model");
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildF2E(const MeshModel &ARefModel) const
    {
        //we just have F2N and E2N
        std::map<TCellID,std::vector<TCellID> > tmp_N2E;
        for(auto e_id : m_mesh->edges()){
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<TCellID> e_node_ids;
            e.getIDs<Node>(e_node_ids);
            for(auto n_id : e_node_ids) {
                tmp_N2E[n_id].push_back(e_id);
            }
        }

        for(auto f_id : m_mesh->faces()){
            Face f = m_mesh->get<Face>(f_id);
            std::vector<TCellID> f_node_ids;
            f.getIDs<Node>(f_node_ids);
            std::vector<TCellID> multiset_edge;
            for(auto n_id : f_node_ids) {
                std::vector<TCellID> current_edges= tmp_N2E[n_id];
                multiset_edge.insert(multiset_edge.end(),current_edges.begin(),current_edges.end());
            }
            std::vector<TCellID> result = keepFilter(multiset_edge,2);
            f.set<Edge>(result);
        }
    }
    /*----------------------------------------------------------------------------*/
    void MeshDoctor::buildR2E(const MeshModel &ARefModel) const
    {
        //we just have R2N and E2N
        std::map<TCellID,std::vector<TCellID> > tmp_N2E;
        for(auto e_id : m_mesh->edges()){
            Edge e = m_mesh->get<Edge>(e_id);
            std::vector<TCellID> e_node_ids;
            e.getIDs<Node>(e_node_ids);
            for(auto n_id : e_node_ids) {
                tmp_N2E[n_id].push_back(e_id);
            }
        }

        for(auto r_id : m_mesh->regions()){
            Region r = m_mesh->get<Region>(r_id);
            std::vector<TCellID> r_node_ids;
            r.getIDs<Node>(r_node_ids);
            std::vector<TCellID> multiset_edges;
            for(auto n_id : r_node_ids) {
                std::vector<TCellID> current_edges = tmp_N2E[n_id];
                multiset_edges.insert(multiset_edges.end(),current_edges.begin(),current_edges.end());
            }
            std::vector<TCellID> r_face_ids = keepFilter(multiset_edges,2);

            r.set<Edge>(r_face_ids);

        }

    }
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
