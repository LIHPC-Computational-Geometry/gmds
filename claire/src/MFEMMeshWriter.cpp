//
// Created by rochec on 26/09/2023.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/MFEMMeshWriter.h>
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
	m_manager->initAndLinkFrom2DMesh(m_mesh, m_linker);

	// Build edges and connectivities
	MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();	// Orient all the faces. WARNING: The orientation is important in .mesh file for MFEM.
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
	writeElements();
	writeBnd();
	writeNodes();

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
	m_stream << "MFEM NC mesh v1.0\n";
	//m_stream << "MFEM mesh v1.0\n";
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
	m_stream << "2\n";
	m_stream << "\n";
}
/*------------------------------------------------------------------------*/
void MFEMMeshWriter::writeElements()
{
	m_stream << "# rank attr geom ref_type nodes/children\n";
	m_stream << "elements\n";
	m_stream << m_mesh->getNbFaces() << "\n";

	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();
		int geom(0);
		int attrib(0);
		if (f_nodes.size()==3){geom=2;attrib=1;}			// geom= TRIANGLE
		else if (f_nodes.size()==4){geom=3;attrib=2;}	// geom= QUAD
		m_stream << 0 << " ";				// RANK
		m_stream << attrib << " ";		// Attribute
		m_stream << geom << " ";	// Geom
		m_stream << 0;		// Ref type
		if (geom==2)
		{
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else if (geom==3)
		{
			m_stream << " " << m_ordering_nodes[f_nodes[0].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[3].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[2].id()];
			m_stream << " " << m_ordering_nodes[f_nodes[1].id()];
		}
		else
		{
			for (auto const n:f_nodes)
			{
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
		/*
		// Test with control on the boundary labeling
		int geom_id(0);
		if (abs(e_nodes[0].X() + 100) < pow(10,-6)
		    && abs(e_nodes[1].X() + 100) < pow(10,-6))
		{
			geom_id = 1;
		}
		else if (abs(e_nodes[0].X() - 200) < pow(10,-6)
		    && abs(e_nodes[1].X() - 200) < pow(10,-6))
		{
			geom_id = 2;
		}
		else if (abs(e_nodes[0].Y() - 170) < pow(10,-6)
		    && abs(e_nodes[1].Y() - 170) < pow(10,-6))
		{
			geom_id = 3;
		}
		else if (abs(e_nodes[0].Y() + 80) < pow(10,-6)
		         && abs(e_nodes[1].Y() + 80) < pow(10,-6))
		{
			geom_id = 4;
		}
		*/
		m_stream << geom_id << " ";	// Attribute
		m_stream << 1 << " ";	// Geom (1= SEGMENT)
		m_stream << m_ordering_nodes[e_nodes[0].id()] << " " << m_ordering_nodes[e_nodes[1].id()] << "\n";
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
