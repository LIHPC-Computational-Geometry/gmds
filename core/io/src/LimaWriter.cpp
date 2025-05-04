/*----------------------------------------------------------------------------*/
#include "gmds/io/LimaWriter.h"
#include "gmds/io/LimaWriterStreamMli2.h"
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <Lima/lima++.h>
#include <LimaP/reader.h>
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
LimaWriter::LimaWriter(gmds::Mesh& AMesh)
  :m_mesh(AMesh), m_nodes_connection(1), m_length_unit(1.)
{}
/*----------------------------------------------------------------------------*/
LimaWriter::~LimaWriter()
{
}
/*----------------------------------------------------------------------------*/
void LimaWriter::setLengthUnit(double AUnit)
{
	m_length_unit = AUnit;
}
/*----------------------------------------------------------------------------*/
void LimaWriter::write(const std::string& AFileName, gmds::MeshModel AModel, int ACompact)
{
	/* Detection du format par le suffixe du nom du fichier. */
	Lima::format_t format = Lima::SUFFIXE;
	format = Lima::_Reader::detectFormat(AFileName);

	switch(format) {
	case Lima::SUFFIXE :
		throw GMDSException("GMDSCEAWriter::write Extension de fichier inconnue\n");
		break;
	case Lima::INCONNU :
		throw GMDSException("GMDSCEAWriter::write Extension de fichier inconnue\n");
		break;
	case Lima::MALIPP2 :
	{
		try {
			gmds::LimaWriterStreamMli2 w(m_mesh);
			w.setLengthUnit(m_length_unit);
			w.write(AFileName,AModel,ACompact);
		}
		catch(gmds::GMDSException& e) {
			std::cerr<<"LimaWriter ERREUR : "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}
	break;
	default :
	{
		try {
			this->writeClassic(AFileName,AModel,ACompact);
		}
		catch(gmds::GMDSException& e) {
			std::cerr<<"LimaWriter ERREUR : "<<e.what()<<std::endl;
			throw GMDSException(e.what());
		}
	}
	break;
	}
}
/*----------------------------------------------------------------------------*/
void LimaWriter::writeClassic(const std::string& AFileName, gmds::MeshModel AModel, int ACompact)
{
	m_nodes_connection.clear();

	TCellID max_id = 0;
	for(auto i: m_mesh.nodes())
	{
		gmds::Node n = m_mesh.get<gmds::Node>(i);
		if(n.id()>max_id)
			max_id = n.id();
	}
	m_nodes_connection.resize(max_id+1);

	Lima::Maillage m;
	m.unite_longueur(m_length_unit);

	writeNodes(m);

	if (AModel.has(E) && m_mesh.getModel().has(E))
		writeEdges(m);

	if (AModel.has(F) && m_mesh.getModel().has(F))
		writeFaces(m);

	if (AModel.has(R) && m_mesh.getModel().has(R))
		writeRegions(m);

	try{
		if(m_mesh.getModel().has(DIM2)){
			m.dimension(Lima::D2);
		}
		m.ecrire(AFileName, Lima::SUFFIXE,1,ACompact);
	}
	catch(Lima::write_erreur& e)
	{
		std::cerr<<"LIMA ERREUR : "<<e.what()
		          <<"\nCela peut venir d'un chemin incorrect, d'un problème de permissions ou de quota."
		          <<std::endl;
		throw GMDSException(e.what());
	}
}
/*----------------------------------------------------------------------------*/
void LimaWriter::writeNodes(Lima::Maillage& ALimaMesh)
{

	for(auto i: m_mesh.nodes())
	{
		gmds::Node n = m_mesh.get<gmds::Node>(i);
		Lima::Noeud n2(n.id()+1,n.X(), n.Y(), n.Z());
		/* we keep the connection between lima node and gmds node through local ids.
		 */
		ALimaMesh.ajouter(n2);
		m_nodes_connection[n.id()] = n2;

	}

	for(auto it = m_mesh.groups_begin<Node>(); it!=m_mesh.groups_end<Node>(); it++)
	{
		Lima::Nuage lima_cl((*it)->name());
		ALimaMesh.ajouter(lima_cl);

		auto ids = (*it)->cells();
		for(unsigned int index =0; index <ids.size(); index++)
			lima_cl.ajouter(m_nodes_connection[ids[index]]);
	}
}
/*----------------------------------------------------------------------------*/
void LimaWriter::writeEdges(Lima::Maillage& ALimaMesh)
{

	for(auto i: m_mesh.edges())
	{
		gmds::Edge e = m_mesh.get<gmds::Edge>(i);
		std::vector<TCellID> nodes = e.getIDs<Node>();
		Lima::Noeud n1 = m_nodes_connection[nodes[0]];
		Lima::Noeud n2 = m_nodes_connection[nodes[1]];
		Lima::Bras e2(e.id()+1,n1, n2);
		ALimaMesh.ajouter(e2);
	}

	for(auto it = m_mesh.groups_begin<Edge>(); it!=m_mesh.groups_end<Edge>(); it++)
	{
		Lima::Ligne lima_li((*it)->name());
		ALimaMesh.ajouter(lima_li);

		auto ids = (*it)->cells();
		for(unsigned int index=0; index<ids.size();index++){
			Lima::Bras b = ALimaMesh.bras_id(ids[index]+1);
			lima_li.ajouter(b);
		}
	}
}
/*----------------------------------------------------------------------------*/
void LimaWriter::writeFaces(Lima::Maillage& ALimaMesh)
{
	for(auto i: m_mesh.faces())
	{
		gmds::Face f = m_mesh.get<gmds::Face>(i);
		std::vector<TCellID> nodes = f.getIDs<Node>();
		switch(f.type()){
		case GMDS_QUAD:{
			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];
			Lima::Polygone f2(f.id()+1,n1, n2,n3,n4);
			ALimaMesh.ajouter(f2);}
		break;
		case GMDS_TRIANGLE:{
			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];

			Lima::Polygone f2(f.id()+1,n1, n2,n3);
			ALimaMesh.ajouter(f2);}
		break;
		case GMDS_POLYGON:{
			switch(nodes.size()){
			case 3:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Polygone f2(f.id()+1,n1, n2, n3);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 4:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Polygone f2(f.id()+1,n1, n2, n3, n4);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 5:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4]];
				Lima::Polygone f2(f.id()+1,n1, n2, n3, n4, n5);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 6:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4]];
				Lima::Noeud n6 = m_nodes_connection[nodes[5]];
				Lima::Polygone f2(f.id()+1,n1, n2, n3, n4, n5, n6);
				ALimaMesh.ajouter(f2);
			}
			break;
			default:
				std::cout<<"Unable to convert a polygon with more than 6 nodes in Lima format"<<std::endl;
			};
		}
		break;
		default:
			std::cout<<"Unable to convert this type of cell"<<std::endl;
		};

	}

	for(auto it = m_mesh.groups_begin<Face>(); it!=m_mesh.groups_end<Face>(); it++)
	{
		Lima::Surface lima_surf((*it)->name());
		ALimaMesh.ajouter(lima_surf);

		auto ids = (*it)->cells();
		for(unsigned int face_index=0; face_index<ids.size();face_index++)
		{
			Lima::Polygone p = ALimaMesh.polygone_id(ids[face_index]+1);
			lima_surf.ajouter(p);
		}
	}


}
/*----------------------------------------------------------------------------*/
void LimaWriter::writeRegions(Lima::Maillage& ALimaMesh)
{
	for(auto i: m_mesh.regions())
	{
		gmds::Region r = m_mesh.get<gmds::Region>(i);
		std::vector<TCellID> nodes = r.getIDs<Node>();
		switch(r.type()){
		case GMDS_TETRA:{
			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];

			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4);
			ALimaMesh.ajouter(f2);
		}
		break;
		case GMDS_HEX:{

			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];
			Lima::Noeud n5 = m_nodes_connection[nodes[4]];
			Lima::Noeud n6 = m_nodes_connection[nodes[5]];
			Lima::Noeud n7 = m_nodes_connection[nodes[6]];
			Lima::Noeud n8 = m_nodes_connection[nodes[7]];


			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4,n5,n6,n7,n8);
			ALimaMesh.ajouter(f2);
		}
		break;
		case GMDS_PYRAMID:{

			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];
			Lima::Noeud n5 = m_nodes_connection[nodes[4]];


			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4,n5);
			ALimaMesh.ajouter(f2);
		}
		break;
		case GMDS_PRISM3:{

			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];
			Lima::Noeud n5 = m_nodes_connection[nodes[4]];
			Lima::Noeud n6 = m_nodes_connection[nodes[5]];

			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4,n5,n6);
			ALimaMesh.ajouter(f2);
		}
		break;

		case GMDS_PRISM5:{

			Lima::Noeud n1 = m_nodes_connection[nodes[0]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3]];
			Lima::Noeud n5 = m_nodes_connection[nodes[4]];
			Lima::Noeud n6 = m_nodes_connection[nodes[5]];
			Lima::Noeud n7 = m_nodes_connection[nodes[6]];
			Lima::Noeud n8 = m_nodes_connection[nodes[7]];
			Lima::Noeud n9 = m_nodes_connection[nodes[8]];
			Lima::Noeud n10= m_nodes_connection[nodes[9]];

			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10);
			ALimaMesh.ajouter(f2);
		}
		break;

		case GMDS_PRISM6:{

			Lima::Noeud n1 = m_nodes_connection[nodes[0 ]];
			Lima::Noeud n2 = m_nodes_connection[nodes[1 ]];
			Lima::Noeud n3 = m_nodes_connection[nodes[2 ]];
			Lima::Noeud n4 = m_nodes_connection[nodes[3 ]];
			Lima::Noeud n5 = m_nodes_connection[nodes[4 ]];
			Lima::Noeud n6 = m_nodes_connection[nodes[5 ]];
			Lima::Noeud n7 = m_nodes_connection[nodes[6 ]];
			Lima::Noeud n8 = m_nodes_connection[nodes[7 ]];
			Lima::Noeud n9 = m_nodes_connection[nodes[8 ]];
			Lima::Noeud n10= m_nodes_connection[nodes[9 ]];
			Lima::Noeud n11= m_nodes_connection[nodes[10]];
			Lima::Noeud n12= m_nodes_connection[nodes[11]];

			Lima::Polyedre f2(r.id()+1,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12);
			ALimaMesh.ajouter(f2);
		}
		break;
		case GMDS_POLYHEDRA:{
			switch(nodes.size()){
			case 4:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 5:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4, n5);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 6:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4]];
				Lima::Noeud n6 = m_nodes_connection[nodes[5]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4, n5, n6);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 8:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4]];
				Lima::Noeud n6 = m_nodes_connection[nodes[5]];
				Lima::Noeud n7 = m_nodes_connection[nodes[6]];
				Lima::Noeud n8 = m_nodes_connection[nodes[7]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4 , n5 , n6, n7, n8);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 10:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0 ]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1 ]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2 ]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3 ]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4 ]];
				Lima::Noeud n6 = m_nodes_connection[nodes[5 ]];
				Lima::Noeud n7 = m_nodes_connection[nodes[6 ]];
				Lima::Noeud n8 = m_nodes_connection[nodes[7 ]];
				Lima::Noeud n9 = m_nodes_connection[nodes[8 ]];
				Lima::Noeud n10= m_nodes_connection[nodes[9 ]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4 , n5 , n6, n7, n8, n9, n10);
				ALimaMesh.ajouter(f2);
			}
			break;
			case 12:
			{
				Lima::Noeud n1 = m_nodes_connection[nodes[0 ]];
				Lima::Noeud n2 = m_nodes_connection[nodes[1 ]];
				Lima::Noeud n3 = m_nodes_connection[nodes[2 ]];
				Lima::Noeud n4 = m_nodes_connection[nodes[3 ]];
				Lima::Noeud n5 = m_nodes_connection[nodes[4 ]];
				Lima::Noeud n6 = m_nodes_connection[nodes[5 ]];
				Lima::Noeud n7 = m_nodes_connection[nodes[6 ]];
				Lima::Noeud n8 = m_nodes_connection[nodes[7 ]];
				Lima::Noeud n9 = m_nodes_connection[nodes[8 ]];
				Lima::Noeud n10= m_nodes_connection[nodes[9 ]];
				Lima::Noeud n11= m_nodes_connection[nodes[10]];
				Lima::Noeud n12= m_nodes_connection[nodes[11]];
				Lima::Polyedre f2(r.id()+1,n1, n2, n3, n4 , n5 , n6,
				                  n7, n8, n9, n10, n11, n12);
				ALimaMesh.ajouter(f2);
			}
			break;
			default:
				std::cout<<"Unable to convert a polyhedron with 7, 9 or more than 12 nodes in Lima format"<<std::endl;
			}
		}
		break;
		default:
			std::cout<<"Unable to convert a polyhedron with 7, 9 or more than 12 nodes in Lima format"<<std::endl;
		}
	}

	for(auto it = m_mesh.groups_begin<Region>(); it!=m_mesh.groups_end<Region>(); it++)
	{
		Lima::Volume lima_vol((*it)->name());
		ALimaMesh.ajouter(lima_vol);

		auto ids = (*it)->cells();
		for(unsigned int region_index=0; region_index<ids.size();region_index++)
		{
			Lima::Polyedre p = ALimaMesh.polyedre_id(ids[region_index]+1);
			lima_vol.ajouter(p);
		}
	}
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/