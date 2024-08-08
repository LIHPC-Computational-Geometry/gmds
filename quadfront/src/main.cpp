#include <iostream>
#include <gmds/quadfront/Quadfront.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>

/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
namespace quadfront {

void main()
{
	std::cout << "============== Début Maillage par avancée de front ==============" << std::endl;

	quadfront::Quadfront mesh("/home/pagea/Documents/Travail/data/cercle2.vtk");

	Variable<int> *front = mesh.m_mesh->newVariable<int, GMDS_EDGE>("front");

	Variable<int> *quad_var = mesh.m_mesh->newVariable<int, GMDS_FACE>("quad");

	std::cout << "liste : " << mesh.m_listBoundary[0].size() << std::endl;
	std::cout << mesh.m_nodeBoundary.size() << std::endl;

	int pr = 10;
	int p = 1;

	while (pr <= 50) {
		std::vector<std::vector<Edge>> listBoundary;
		for (auto &l : mesh.m_listBoundary) {
			std::vector<Edge> listEdge;
			for (auto &e : l) {
				if (mesh.m_nodeBoundary.find(e.get<Node>()[0]) != mesh.m_nodeBoundary.end()
				    and mesh.m_nodeBoundary.find(e.get<Node>()[1]) != mesh.m_nodeBoundary.end()) {
					Face quad = mesh.edgeRecovery(e);
					std::map<Node, std::pair<Edge, Edge>> mapQuad = mesh.get_edgeInFace(quad);

					Edge edge_up;
					if (std::get<0>(mapQuad[e.get<Node>()[0]]) == e) {
						Node node_opposit0 = std::get<1>(mapQuad[e.get<Node>()[0]]).getOppositeNode(e.get<Node>()[0]);
						edge_up = std::get<0>(mapQuad[node_opposit0]) == std::get<1>(mapQuad[e.get<Node>()[0]]) ? std::get<1>(mapQuad[node_opposit0])
						                                                                                        : std::get<0>(mapQuad[node_opposit0]);
					}
					else {
						Node node_opposit0 = std::get<0>(mapQuad[e.get<Node>()[0]]).getOppositeNode(e.get<Node>()[0]);
						edge_up = std::get<0>(mapQuad[node_opposit0]) == std::get<0>(mapQuad[e.get<Node>()[0]]) ? std::get<1>(mapQuad[node_opposit0])
						                                                                                        : std::get<0>(mapQuad[node_opposit0]);
					}
					gmds::IGMeshIOService ioService(mesh.m_mesh);
					gmds::VTKWriter vtkWriter(&ioService);
					vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
					vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
					vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link2.vtk");

					listEdge.push_back(edge_up);
					front->set(edge_up.id(), pr);
					quad_var->set(quad.id(), p);

					mesh.interiorSmoothing(edge_up.get<Node>()[0]);
					mesh.interiorSmoothing(edge_up.get<Node>()[1]);
					if (mesh.m_edgeBoundary[e] > mesh.m_listBoundary.size()) {
						mesh.interiorSmoothing(e.get<Node>()[0]);
						mesh.interiorSmoothing(e.get<Node>()[1]);
					}
				}
			}
			for (auto &ex : listEdge) {
				std::cout << " front " << ex.id() << " : " << ex.get<Node>()[0] << " ; " << ex.get<Node>()[1] << std::endl;
			}
			for (auto &no : mesh.m_nodeBoundary) {
				std::cout << " Node " << no.first << std::endl;
				std::cout << "\t Edge " << std::get<1>(no.second).id() << " : " << std::get<1>(no.second).get<Node>()[0] << " ; "
				          << std::get<1>(no.second).get<Node>()[1] << std::endl;
				std::cout << "\t Edge " << std::get<2>(no.second).id() << " : " << std::get<2>(no.second).get<Node>()[0] << " ; "
				          << std::get<2>(no.second).get<Node>()[1] << std::endl;
			}
			listBoundary.push_back(listEdge);
		}

		mesh.m_listBoundary = listBoundary;
		pr += 10;
		p += 1;

		for (auto &m : mesh.m_listBoundary) {
			for (auto &ef : m) {
				std::cout << "\n=== SMOOTH EDGE " << ef.id() << "===\n " << std::endl;
				Face quad = ef.get<Face>()[0].nbNodes() == 4 ? ef.get<Face>()[0] : ef.get<Face>()[1];
				Edge edge_base = mesh.get_oppositEdge_Quad(ef, quad);
				std::cout << "ef : " << ef.get<Node>()[0] << " ; " << ef.get<Node>()[1] << std::endl;
				std::cout << "edge_base : " << edge_base.get<Node>()[0] << " ; " << edge_base.get<Node>()[1] << std::endl;
				if (mesh.m_edgeBoundary[edge_base] > mesh.m_listBoundary.size()) {
					mesh.interiorSmoothing(edge_base.get<Node>()[0]);
					mesh.interiorSmoothing(edge_base.get<Node>()[1]);
				}

				std::cout << "\n=== SMOOTH Node " << ef.get<Node>()[0].id() << " ===\n " << std::endl;
				mesh.frontSmoothing(ef.get<Node>()[0]);
				std::cout << "\n=== SMOOTH Node " << ef.get<Node>()[1].id() << " ===\n " << std::endl;
				mesh.frontSmoothing(ef.get<Node>()[1]);

				/*
				for (auto &node : ef.get<Node>()) {
				   for (auto &e0 : node.get<Edge>()) {
				      if (mesh.m_edgeBoundary.find(e0) == mesh.m_edgeBoundary.end()) {
				         if (e0.get<Face>()[0].nbEdges() == 3 and e0.get<Face>()[1].nbEdges() == 3) {
				            Node opposit_node = e0.getOppositeNode(node);
				            static std::vector<math::Point> vec_point;
				            for (auto &ee0 : opposit_node.get<Edge>()) {
				               Node node_for = ee0.getOppositeNode(opposit_node);
				               vec_point.push_back(node_for.point());
				            }
				            math::Point center = math::Point::massCenter(vec_point);
				            e0.getOppositeNode(node).X() = center.X();
				            e0.getOppositeNode(node).Y() = center.Y();
				            e0.getOppositeNode(node).Y() = center.Y();
				         }
				      }
				   }
				}
				 */
			}
			std::cout << "suivant" << std::endl;
		}

		gmds::IGMeshIOService ioService(mesh.m_mesh);
		gmds::VTKWriter vtkWriter(&ioService);
		vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
		vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link2.vtk");
	}
}


/*----------------------------------------------------------------------------*/
}     // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/