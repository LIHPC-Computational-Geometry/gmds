//
// Created by rochec on 26/09/2023.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/MFEMMeshWriter.h>
#include <gmds/claire/Utils.h>
#include <gmds/ig/MeshDoctor.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
MFEMMeshWriter::MFEMMeshWriter(Mesh *AMesh, std::string AFileName) :
	m_mesh(AMesh),
  m_filename(AFileName),
  m_manager(new cad::FACManager()),
  m_linker(new cad::GeomMeshLinker())
{
	m_stream= std::ofstream(AFileName+".mesh", std::ios::out);

	// Build edges and connectivities
	if (m_mesh->getNbRegions()==0)
	{
		m_manager->initAndLinkFrom2DMesh(m_mesh, m_linker);
		MeshDoctor doc(m_mesh);
		doc.buildEdgesAndX2E();
		doc.updateUpwardConnectivity();
		doc.orient2DFaces();     // Orient all the faces. WARNING: The orientation is important in .mesh file for MFEM.
	}
	else
	{
		MeshDoctor doc(m_mesh);
		doc.buildFacesAndR2F();
		doc.updateUpwardConnectivity();
		doc.orient2DFaces();

		for (auto r_id:m_mesh->regions())
		{
			math::Utils::orientRegion(m_mesh, m_mesh->get<Region>(r_id));
		}
	}
}


/*------------------------------------------------------------------------*/
MFEMMeshWriter::STATUS MFEMMeshWriter::execute()
{
	std::cout << "MFEM WRITER ONLY FOR 2D..." << std::endl;	// Because the dimension is set to 2 in header() method

	// Build edges and connectivities
	/*
	gmds::MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	 */

	// Build the ordering node map
	orderingNodes();

	// .mesh writing for MFEM
	m_stream.precision(15);

	header();
	if (m_mesh->getNbRegions()==0) 	// 3D writer
	{
		std::cout << "2D writing..." << std::endl;
		writeElements();
		writeBnd();
		writeNodes();
	}
	else	// 2D writer
	{
		std::cout << "3D writing..." << std::endl;
		writeElements3D();
		writeBnd3D();
		writeNodes3D();
	}

	m_stream << "mfem_mesh_end\n";

	return MFEMMeshWriter::SUCCESS;
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::orderingNodes()
{
	int compteur(0);
	for (auto n_id:m_mesh->nodes())
	{
		m_ordering_nodes[n_id] = compteur;
		compteur++;
	}
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::header()
{
	// Header
	//m_stream << "MFEM NC mesh v1.0\n";
	m_stream << "MFEM mesh v1.0\n";
	m_stream << "\n";
	m_stream << "# NCMesh supported geometry types:\n";
	m_stream << "# SEGMENT     = 1\n";
	m_stream << "# TRIANGLE    = 2\n";
	m_stream << "# SQUARE      = 3\n";
	m_stream << "# TETRAHEDRON = 4\n";
	m_stream << "# CUBE        = 5\n";
	m_stream << "# PRISM       = 6\n";
	m_stream << "# PYRAMID     = 7\n";
	m_stream << "\n";
	m_stream << "dimension\n";
	if (m_mesh->getNbRegions()==0)
	{
		m_stream << "2\n";
	}
	else
	{
		m_stream << "3\n";
	}
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeElements()
{
	m_stream << "# rank attr geom ref_type nodes/children\n";
	m_stream << "elements\n";
	m_stream << m_mesh->getNbFaces() << "\n";

	for (auto f_id : m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		int geom(0);
		int attrib(0);
		if (f_nodes.size() == 3) {
			geom = 2;
			attrib = 1;
		}     // geom= TRIANGLE
		else if (f_nodes.size() == 4) {
			geom = 3;
			attrib = 2;
		}                              // geom= QUAD
		m_stream << 0 << " ";          // RANK
		m_stream << attrib << " ";     // Attribute
		m_stream << geom << " ";       // Geom
		m_stream << 0;                 // Ref type
		if (geom == 2) {
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else if (geom == 3) {
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[3].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else {
			for (auto const n : f_nodes) {
				m_stream << " " << m_ordering_nodes[n.id()];
			}
		}

		m_stream << "\n";
	}
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeElements3D()
{
	m_stream << "# rank attr geom ref_type nodes/children\n";
	m_stream << "elements\n";
	m_stream << m_mesh->getNbRegions() << "\n";
	for (auto r_id:m_mesh->regions())
	{
		Region r = m_mesh->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		int geom(0);
		int attrib(0);
		if (r.type() == GMDS_HEX) // geom = HEXA
		{
			geom = 5;
			attrib = 1;
		}
		else if (r.type() == GMDS_PYRAMID) // geom = PYRAMID
		{
			geom = 7;
			attrib = 2;
		}
		else if (r.type() == GMDS_TETRA) // geom = TETRA
		{
			geom = 4;
			attrib = 1;
		}
		else
		{
			std::cout << "WARNING MFEMMeshWriter: 3D Region not implemented yet." << std::endl;
		}
		//m_stream << 0 << " ";          // RANK
		m_stream << attrib << " ";     // Attribute
		m_stream << geom << " ";       // Geom
		//m_stream << 0;
		if (geom == 5)
		{
			m_stream << " " << m_ordering_nodes[r_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[3].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[1].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[4].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[7].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[6].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[5].id()];
		}
		else if (geom == 7)
		{
			m_stream << " " << m_ordering_nodes[r_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[3].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[1].id()];
			m_stream << " " << m_ordering_nodes[r_nodes[4].id()];
		}
		else {
			for (auto const n : r_nodes) {
				m_stream << " " << m_ordering_nodes[n.id()];
			}
		}
		m_stream << "\n";
	}
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeBnd()
{
	std::vector<Edge> bnd_edges;
	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		if (e.get<Face>().size() == 1)
		{
			bnd_edges.push_back(e);
		}
	}

	m_stream << "# attr geom nodes\n";
	m_stream << "boundary\n";
	m_stream << bnd_edges.size() << "\n";
	for (auto e:bnd_edges)
	{
		std::vector<Node> e_nodes = e.get<Node>();
		int geom_id = m_linker->getGeomId(e)+1;
		m_stream << geom_id << " ";	// Attribute
		m_stream << 1 << " ";	// Geom (1= SEGMENT)
		m_stream << m_ordering_nodes[e_nodes[0].id()] << " " << m_ordering_nodes[e_nodes[1].id()] << "\n";
	}
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeBnd3D()
{
	std::vector<Face> bnd_faces;
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		if (f.get<Region>().size() == 1)
		{
			bnd_faces.push_back(f);
		}
	}

	m_stream << "# attr geom nodes\n";
	m_stream << "boundary\n";
	m_stream << bnd_faces.size() << "\n";
	std::cout << "Nbr bnd faces: " << bnd_faces.size() << std::endl;
	for (auto f:bnd_faces)
	{
		std::vector<Node> f_nodes = f.get<Node>();
		int geom(0);
		int attrib(0);
		if (f_nodes.size() == 3)
		{
			geom = 2;
			attrib = 1;
		}     // geom= TRIANGLE
		else if (f_nodes.size() == 4)
		{
			geom = 3;
			attrib = 2;
		}                              // geom= QUAD
		m_stream << attrib << " ";     // Attribute
		m_stream << geom ;       // Geom
		if (geom == 2) {
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else if (geom == 3) {
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[3].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else {
			for (auto const n : f_nodes) {
				m_stream << " " << m_ordering_nodes[n.id()];
			}
		}
		m_stream << "\n";
	}
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeNodes()
{
	/*
	m_stream << "vertices\n";
	m_stream << m_mesh->getNbNodes() << "\n";
	m_stream << "2\n";
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_stream << n.X() << " " << n.Y() << "\n";
	}
	*/

	m_stream << "# mesh curvature GridFunction\n";
	m_stream << "nodes\n";
	m_stream << "FiniteElementSpace\n";
	m_stream << "FiniteElementCollection: H1_2D_P1\n";
	m_stream << "VDim: 2\n";
	m_stream << "Ordering: 1\n";
	m_stream << "\n";
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_stream << n.X() << " " << n.Y() << "\n";
	}

	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeNodes3D()
{
	m_stream << "vertices\n";
	m_stream << m_mesh->getNbNodes() << "\n";
	m_stream << "\n";
	//m_stream << "# mesh curvature GridFunction\n";
	m_stream << "nodes\n";
	m_stream << "FiniteElementSpace\n";
	m_stream << "FiniteElementCollection: H1_3D_P1\n";
	m_stream << "VDim: 3\n";
	m_stream << "Ordering: 1\n";
	m_stream << "\n";
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_stream << n.X() << " " << n.Y() << " " << n.Z() << "\n";
	}

	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
