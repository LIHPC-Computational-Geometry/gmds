//
// Created by rochec on 11/01/24.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/SU2Writer_3D.h>
#include <gmds/aero/AeroBoundaries_3D.h>
#include <gmds/aero/Utils.h>
#include <gmds/ig/MeshDoctor.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
SU2Writer_3D::SU2Writer_3D(Mesh *AMesh, std::string AFileName) :
	m_mesh(AMesh),
	m_filename(AFileName)
{
	MeshDoctor doc(m_mesh);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}


/*------------------------------------------------------------------------*/
SU2Writer_3D::STATUS SU2Writer_3D::execute()
{

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream(m_filename, std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);

	math::Utils::MeshCleaner(m_mesh);

	std::map<TCellID,TCellID> map_new_id;
	int iterateur = 0;
	for (auto n_id:m_mesh->nodes())
	{
		map_new_id[n_id] = iterateur ;
		iterateur++;
	}

	//Problem dimension
	stream << "%\n";
	stream << "% Problem dimension \n";
	stream << "%\n";
	stream << "NDIME= " << 3 << "\n";

	// Liste des éléments
	stream << "%\n";
	stream << "% Inner element connectivity \n";
	stream << "%\n";
	stream << "NELEM= " << m_mesh->getNbRegions() << "\n";
	for (auto r_id:m_mesh->regions())
	{
		Region r = m_mesh->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		if (r_nodes.size() == 4){	// TETRA
			stream << 10 << " " << map_new_id[r_nodes[0].id()] << " " << map_new_id[r_nodes[1].id()] <<
			   " " << map_new_id[r_nodes[2].id()] << " " << map_new_id[r_nodes[3].id()] << "\n";
		}
		else if (r_nodes.size() == 8){	// HEX
			stream << 12 << " " << map_new_id[r_nodes[0].id()] <<
			   " " << map_new_id[r_nodes[1].id()] <<
			   " " << map_new_id[r_nodes[2].id()] <<
			   " " << map_new_id[r_nodes[3].id()] <<
			   " " << map_new_id[r_nodes[4].id()] <<
			   " " << map_new_id[r_nodes[5].id()] <<
			   " " << map_new_id[r_nodes[6].id()] <<
			   " " << map_new_id[r_nodes[7].id()] <<
			   "\n";
		}
		else
		{
			std::cout << "ATTENTION SU2Writer_3D: Element not implemented yet." << std::endl;
			std::cout << "Info: Element has " << r_nodes.size() << " nodes." << std::endl;
		}
	}


	// Liste des noeuds
	stream << "%\n";
	stream << "% Node coordinates \n";
	stream << "%\n";
	stream << "NPOIN= " << m_mesh->getNbNodes() << "\n";

	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		stream << n.X() << " " << n.Y() << " " << n.Z() << " " << map_new_id[n_id] << "\n";

	}



	// Label des bords
	AeroBoundaries_3D* bnd_3D = new AeroBoundaries_3D(m_mesh);
	bnd_3D->execute();

	stream << "%\n";
	stream << "% Boundary elements \n";
	stream << "%\n";
	stream << "NMARK= " << bnd_3D->getNbrBords() << "\n";

	// FARFIELD
	int color_ext = bnd_3D->getColorAmont();
	std::vector<TCellID> faces_id_ff = bnd_3D->BndFaces(color_ext);
	stream << "MARKER_TAG= FARFIELD\n";
	stream << "MARKER_ELEMS= " << faces_id_ff.size() << "\n";
	for (auto f_id : faces_id_ff)
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<Node> face_nodes = f.get<Node>();
		if (face_nodes.size() == 4)
		{
			stream << 9 << " " << map_new_id[face_nodes[0].id()] << " " << map_new_id[face_nodes[1].id()] << " " << map_new_id[face_nodes[2].id()] << " "
			       << map_new_id[face_nodes[3].id()] << "\n";
		}
		}


	// Autres bords
	int wall = 1;
	for (int i=1;i <= bnd_3D->getNbrBords();i++)
	{
		if (i != color_ext)
		{
			std::vector<TCellID> faces_id = bnd_3D->BndFaces(i);
			stream << "MARKER_TAG= wall_" << wall << "\n";
			stream << "MARKER_ELEMS= " << faces_id.size() << "\n";
			for (auto f_id : faces_id)
			{
				Face f = m_mesh->get<Face>(f_id);
				std::vector<Node> face_nodes = f.get<Node>();
				if (face_nodes.size() == 4)
				{
					stream << 9 << " " << map_new_id[face_nodes[0].id()] << " " << map_new_id[face_nodes[1].id()] << " " << map_new_id[face_nodes[2].id()] << " "
					       << map_new_id[face_nodes[3].id()] << "\n";
				}
			}
			wall ++;
		}
	}

	stream.close();

	return SU2Writer_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/