/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/io/LimaReader.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
LimaReader::LimaReader(Mesh& AMesh)
  :m_mesh(AMesh), m_lengthUnit(1.)
{}
/*----------------------------------------------------------------------------*/
LimaReader::~LimaReader()= default;
/*----------------------------------------------------------------------------*/
double LimaReader::getLengthUnit(){
	return m_lengthUnit;
}
/*----------------------------------------------------------------------------*/
void LimaReader::read(const std::string& AFileName, MeshModel AModel){
	Lima::Maillage m;
	try{
		m.lire(AFileName);
	} catch(...)
	{
		throw GMDSException("Lima cannot read the file "+AFileName);
	}
	m_lengthUnit = m.unite_longueur();
	readNodes(m);
	if (m_mesh.getModel().has(E) && AModel.has(E))
		readEdges(m);

	if (m_mesh.getModel().has(F) && AModel.has(F))
		readFaces(m);

	if (m_mesh.getModel().has(R) && AModel.has(R))
		readRegions(m);
}
/*----------------------------------------------------------------------------*/
void LimaReader::readNodes(Lima::Maillage& ALimaMesh){
	/** look for the highest node id*/
	size_t max_id=0;
	for(auto i = 0; i<ALimaMesh.nb_noeuds(); i++)
		if(ALimaMesh.noeud(i).id()>max_id)
			max_id=ALimaMesh.noeud(i).id();

	m_lima2gmds_node_ids.resize(max_id);

	for(auto i = 0; i < ALimaMesh.nb_noeuds(); i++){
		Lima::Noeud ni = ALimaMesh.noeud(i);
		Node n = m_mesh.newNode(ni.x(),ni.y(),ni.z());
		m_lima2gmds_node_ids[ni.id()-1] = n.id();
	}

	for(auto index=0;index<ALimaMesh.nb_nuages();index++){
		Lima::Nuage lima_nuage = ALimaMesh.nuage(index);
		int nbNodesInCloud = lima_nuage.nb_noeuds();

		CellGroup<Node>* cl = m_mesh.newGroup<Node>(lima_nuage.nom());
		for(auto node_index = 0; node_index<nbNodesInCloud;node_index++){
			cl->add(m_lima2gmds_node_ids[lima_nuage.noeud(node_index).id()-1]);
		}
	}

}
/*----------------------------------------------------------------------------*/
void LimaReader::readEdges(Lima::Maillage& ALimaMesh){
	int max_id=0;
	for(auto i = 0; i < ALimaMesh.nb_bras(); i++){
		if(ALimaMesh.bras(i).id()>max_id)
			max_id=ALimaMesh.bras(i).id();
	}

	std::vector<Edge> edges_connection;
	edges_connection.resize(max_id);
	for(auto i = 0; i < ALimaMesh.nb_bras(); i++){
		Lima::Bras  b = ALimaMesh.bras(i);
		Edge e = m_mesh.newEdge(m_lima2gmds_node_ids[b.noeud(0).id()-1],
		                        m_lima2gmds_node_ids[b.noeud(1).id()-1]);
		edges_connection[b.id()-1]=e;
	}

	for(auto index=0;index<ALimaMesh.nb_lignes();index++){
		Lima::Ligne lima_ligne = ALimaMesh.ligne(index);
		int nbEdgesInLine= lima_ligne.nb_bras();

		CellGroup<Edge>* li = m_mesh.newGroup<Edge>(lima_ligne.nom());

		for(auto edge_index = 0; edge_index<nbEdgesInLine;edge_index++){
			li->add(edges_connection[lima_ligne.bras(edge_index).id()-1]);
		}
	}
}
/*----------------------------------------------------------------------------*/
void LimaReader::readFaces(Lima::Maillage& ALimaMesh){
	int max_id=0;
	for(auto i = 0; i < ALimaMesh.nb_polygones(); i++){
		if(ALimaMesh.polygone(i).id()>max_id)
			max_id=ALimaMesh.polygone(i).id();
	}

	std::vector<Face> faces_connection;
	faces_connection.resize(max_id);

	for(auto i = 0; i < ALimaMesh.nb_polygones(); i++){
		Lima::Polygone  p = ALimaMesh.polygone(i);
		Face f;
		switch(p.nb_noeuds()){
		case 3:
		{
			f=m_mesh.newTriangle(m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                       m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                       m_lima2gmds_node_ids[p.noeud(2).id()-1]);
		}
		break;
		case 4:
		{
			f=m_mesh.newQuad( m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                   m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                   m_lima2gmds_node_ids[p.noeud(2).id()-1],
			                   m_lima2gmds_node_ids[p.noeud(3).id()-1]);
		}
		break;
		case 5:
		{
			std::vector<TCellID> nodes;
			nodes.resize(5);
			nodes[0] = m_lima2gmds_node_ids[p.noeud(0).id()-1];
			nodes[1] = m_lima2gmds_node_ids[p.noeud(1).id()-1];
			nodes[2] = m_lima2gmds_node_ids[p.noeud(2).id()-1];
			nodes[3] = m_lima2gmds_node_ids[p.noeud(3).id()-1];
			nodes[4] = m_lima2gmds_node_ids[p.noeud(4).id()-1];
			f=m_mesh.newPolygon(nodes);
		}
		break;
		case 6:
		{
			std::vector<TCellID> nodes;
			nodes.resize(6);
			nodes[0] = m_lima2gmds_node_ids[p.noeud(0).id()-1];
			nodes[1] = m_lima2gmds_node_ids[p.noeud(1).id()-1];
			nodes[2] = m_lima2gmds_node_ids[p.noeud(2).id()-1];
			nodes[3] = m_lima2gmds_node_ids[p.noeud(3).id()-1];
			nodes[4] = m_lima2gmds_node_ids[p.noeud(4).id()-1];
			nodes[5] = m_lima2gmds_node_ids[p.noeud(5).id()-1];
			f=m_mesh.newPolygon(nodes);
		}
		break;
		}
		faces_connection[p.id()-1]=f;
	}

	for(auto index=0;index<ALimaMesh.nb_surfaces();index++){
		Lima::Surface lima_surf = ALimaMesh.surface(index);
		int nbFacesInSurf= lima_surf.nb_polygones();
		CellGroup<Face>* su = m_mesh.newGroup<Face>(lima_surf.nom());
		for(auto face_index = 0; face_index<nbFacesInSurf;face_index++){
			su->add(faces_connection[lima_surf.polygone(face_index).id()-1]);
		}
	}
}
/*----------------------------------------------------------------------------*/
void LimaReader::readRegions(Lima::Maillage& ALimaMesh){

	Lima::size_type max_id=0;
	for(auto i = 0; i < ALimaMesh.nb_polyedres(); i++){
		if(ALimaMesh.polyedre(i).id()>max_id)
			max_id=ALimaMesh.polyedre(i).id();
	}

	std::vector<Region> regions_connection;
	regions_connection.resize(max_id);

	for(auto i = 0; i < ALimaMesh.nb_polyedres(); ++i){
		Lima::Polyedre  p = ALimaMesh.polyedre(i);
		Region  r;
		switch(p.nb_noeuds()){
		case 4:
		{
			r=m_mesh.newTet(m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(2).id()-1],
			                 m_lima2gmds_node_ids[p.noeud(3).id()-1]);
		}
		break;
		case 5:
		{
			r=m_mesh.newPyramid(m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                      m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                      m_lima2gmds_node_ids[p.noeud(2).id()-1],
			                     m_lima2gmds_node_ids[p.noeud(3).id()-1],
			                      m_lima2gmds_node_ids[p.noeud(4).id()-1]);
		}
		break;
		case 6:
		{
			r=m_mesh.newPrism3(m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                     m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                     m_lima2gmds_node_ids[p.noeud(2).id()-1],
			                    m_lima2gmds_node_ids[p.noeud(3).id()-1],
			                     m_lima2gmds_node_ids[p.noeud(4).id()-1],
			                     m_lima2gmds_node_ids[p.noeud(5).id()-1]);
		}
		break;
		case 8:
		{
			r=m_mesh.newHex(m_lima2gmds_node_ids[p.noeud(0).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(1).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(2).id()-1],
			                 m_lima2gmds_node_ids[p.noeud(3).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(4).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(5).id()-1],
			                 m_lima2gmds_node_ids[p.noeud(6).id()-1],
			                  m_lima2gmds_node_ids[p.noeud(7).id()-1]);
		}
		break;
		case 10:
		{
			throw GMDSException("Prism5 type not yet implemented");
			//				r=m_mesh.newPrism5(nodes_connection_[p.noeud(0).id()-1],
			//							      nodes_connection_[p.noeud(1).id()-1],
			//							      nodes_connection_[p.noeud(2).id()-1],
			//							      nodes_connection_[p.noeud(3).id()-1],
			//							      nodes_connection_[p.noeud(4).id()-1],
			//							      nodes_connection_[p.noeud(5).id()-1],
			//							      nodes_connection_[p.noeud(6).id()-1],
			//							      nodes_connection_[p.noeud(7).id()-1],
			//							      nodes_connection_[p.noeud(8).id()-1],
			//								  nodes_connection_[p.noeud(9).id()-1]);
		}
		break;
		case 12:
		{
			throw GMDSException("Prism6 type not yet implemented");
			//				r=m_mesh.newPrism6(nodes_connection_[p.noeud(0).id()-1],
			//							      nodes_connection_[p.noeud(1).id()-1],
			//							      nodes_connection_[p.noeud(2).id()-1],
			//							      nodes_connection_[p.noeud(3).id()-1],
			//							      nodes_connection_[p.noeud(4).id()-1],
			//							      nodes_connection_[p.noeud(5).id()-1],
			//							      nodes_connection_[p.noeud(6).id()-1],
			//							      nodes_connection_[p.noeud(7).id()-1],
			//							      nodes_connection_[p.noeud(8).id()-1],
			//							      nodes_connection_[p.noeud(9).id()-1],
			//							      nodes_connection_[p.noeud(10).id()-1],
			//								  nodes_connection_[p.noeud(11).id()-1]);
		}
		break;
		}
		regions_connection[p.id()-1]=r;
	}

	for(auto index=0;index<ALimaMesh.nb_volumes();index++){
		Lima::Volume lima_vol = ALimaMesh.volume(index);
		auto nbRegionsInVol= lima_vol.nb_polyedres();
		CellGroup<Region>* vo = m_mesh.newGroup<Region>(lima_vol.nom());

		for(auto r_index = 0; r_index<nbRegionsInVol;r_index++){
			vo->add(regions_connection[lima_vol.polyedre(r_index).id()-1]);
		}
	}
}