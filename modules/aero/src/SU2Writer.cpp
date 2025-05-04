//
// Created by rochec on 19/05/2022.
//


/*------------------------------------------------------------------------*/
#include <gmds/aero/SU2Writer.h>
#include <gmds/aero/AeroBoundaries_2D.h>
#include <gmds/aero/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
SU2Writer::SU2Writer(Mesh *AMesh, std::string AFileName, double Ax_lim_inout) {
	m_mesh = AMesh;
	m_filename = AFileName;
	m_x_lim_inout = Ax_lim_inout;
}


/*------------------------------------------------------------------------*/
SU2Writer::STATUS SU2Writer::execute()
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
	//stream << "NDIME= " << m_mesh->getDim() << "\n";
	stream << "NDIME= " << 2 << "\n";

	if (m_mesh->getDim()==3){
		std::cout << "ATTENTION :  Ce n'est pas encore implémenté en 3D." << std::endl;
	}

	// Liste des éléments
	stream << "%\n";
	stream << "% Inner element connectivity \n";
	stream << "%\n";
	stream << "NELEM= " << m_mesh->getNbFaces() << "\n";
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<Node> face_nodes = f.get<Node>();
		if (face_nodes.size() == 3){	// TRIANGLES
			stream << 5 << " " << map_new_id[face_nodes[0].id()] << " " << map_new_id[face_nodes[1].id()] << " " << map_new_id[face_nodes[2].id()] << "\n";
		}
		else if (face_nodes.size() == 4){	// QUAD
			stream << 9 << " " << map_new_id[face_nodes[0].id()] << " " << map_new_id[face_nodes[1].id()] << " " << map_new_id[face_nodes[2].id()] << " " << map_new_id[face_nodes[3].id()] << "\n";
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
		stream << n.X() << " " << n.Y() << " " << map_new_id[n_id] << "\n";

	}


	// Label des bords
	AeroBoundaries_2D* bnd_2D = new AeroBoundaries_2D(m_mesh);
	bnd_2D->execute();

	stream << "%\n";
	stream << "% Boundary elements \n";
	stream << "%\n";
	stream << "NMARK= " << bnd_2D->getNbrBords() << "\n";

	// Bord extérieur
	int color_ext = bnd_2D->getColorAmont();
	std::vector<TCellID> edges_id_ext = bnd_2D->BndEdges(color_ext);

	std::vector<TCellID> edges_outlet;
	for (auto e_id:edges_id_ext)
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> edge_nodes = e.get<Node>();
		if (edge_nodes[0].X() > m_x_lim_inout
		    && edge_nodes[1].X() > m_x_lim_inout)
		{
			edges_outlet.push_back(e_id);
		}
	}
	if (edges_outlet.size() != 0) {
		stream << "MARKER_TAG= outlet\n";
		stream << "MARKER_ELEMS= " << edges_outlet.size() << "\n";
		for (auto e_id : edges_outlet) {
			Edge e = m_mesh->get<Edge>(e_id);
			std::vector<Node> edge_nodes = e.get<Node>();
			stream << 3 << " " << map_new_id[edge_nodes[0].id()] << " " << map_new_id[edge_nodes[1].id()] << "\n";
		}
	}

	std::vector<TCellID> edges_inlet;
	for (auto e_id:edges_id_ext)
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> edge_nodes = e.get<Node>();
		if (edge_nodes[0].X() <= m_x_lim_inout
		    || edge_nodes[1].X() <= m_x_lim_inout)
		{
			edges_inlet.push_back(e_id);
		}
	}
	if (edges_inlet.size() != 0) {
		stream << "MARKER_TAG= inlet\n";
		stream << "MARKER_ELEMS= " << edges_inlet.size() << "\n";
		for (auto e_id : edges_inlet) {
			Edge e = m_mesh->get<Edge>(e_id);
			std::vector<Node> edge_nodes = e.get<Node>();
			stream << 3 << " " << map_new_id[edge_nodes[0].id()] << " " << map_new_id[edge_nodes[1].id()] << "\n";
		}
	}




	// Autres bords
	int wall = 1;
	for (int i=1;i <= bnd_2D->getNbrBords();i++) {
		if (i != color_ext) {
			std::vector<TCellID> edges_id = bnd_2D->BndEdges(i);
			stream << "MARKER_TAG= wall_" << wall << "\n";
			stream << "MARKER_ELEMS= " << edges_id.size() << "\n";
			for (auto e_id : edges_id) {
				Edge e = m_mesh->get<Edge>(e_id);
				std::vector<Node> edge_nodes = e.get<Node>();
				stream << 3 << " " << map_new_id[edge_nodes[0].id()] << " " << map_new_id[edge_nodes[1].id()] << "\n";
			}
			wall ++;
		}
	}


	stream.close();

	return SU2Writer::SUCCESS;
}
/*------------------------------------------------------------------------*/