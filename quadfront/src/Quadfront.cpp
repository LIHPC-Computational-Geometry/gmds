/*----------------------------------------------------------------------------*/
#include "gmds/quadfront/Quadfront.h"
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/VTKReader.h"
#include <cmath>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Vector.h>
#include <variant>
#include <vector>

/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace quadfront {
/*----------------------------------------------------------------------------*/
Quadfront::Quadfront(const std::string &AFilename)
{

	m_mesh = new gmds::Mesh(gmds::MeshModel(gmds::DIM2 | gmds::F | gmds::E | gmds::N | gmds::N2E | gmds::N2F | gmds::E2N | gmds::E2F | gmds::F2N | gmds::F2E));

	Variable<int> *nb_boundary = m_mesh->newVariable<int, GMDS_NODE>("nb_boundary");

	std::cout << "============== Start Initialisation ==============" << std::endl;
	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	// vtkReader.setDataOptions(gmds::N | gmds::F);
	vtkReader.read(AFilename);

	gmds::MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	std::cout << "Nombre de noeuds: " << m_mesh->getNbNodes() << std::endl;
	std::cout << "Nombre d'arêtes: " << m_mesh->getNbEdges() << std::endl;
	std::cout << "Nombre de faces: " << m_mesh->getNbFaces() << std::endl;

	std::vector<Edge> FrontEdges;

	for (auto &e : m_mesh->edges()) {
		Edge edge = m_mesh->get<Edge>(e);
		if (edge.nbFaces() == 1) {
			FrontEdges.push_back(edge);
		}
	}

	int i = 1;
	while (!FrontEdges.empty()) {
		std::vector<Edge> listBoundary;
		Edge e0 = FrontEdges.front();
		std::vector<Node> e0_nodes = e0.get<Node>();

		Node nEnd = e0_nodes[0];
		Node nStart = e0_nodes[1];
		m_nodeFront[nEnd] = std::make_tuple(0, e0, e0);
		m_nodeBoundary[nEnd] = std::make_pair(e0, e0);
		listBoundary.push_back(e0);

		while (nStart != nEnd) {
			std::vector<Edge> vEdge = nStart.get<Edge>();

			auto it = std::find(FrontEdges.begin(), FrontEdges.end(), e0);
			if (it != FrontEdges.end()) {
				// Si l'élément est trouvé, le supprimer
				FrontEdges.erase(it);
			}

			std::sort(vEdge.begin(), vEdge.end());
			std::sort(FrontEdges.begin(), FrontEdges.end());

			// Obtenir arête adjacente
			std::vector<Edge> intersection;
			std::set_intersection(vEdge.begin(), vEdge.end(), FrontEdges.begin(), FrontEdges.end(), std::back_inserter(intersection));
			if (size(intersection) == 0) {     // Intersection is empty
				break;
			}
			else {
				m_nodeFront[nStart] = std::make_tuple(0, e0, intersection[0]);
				m_nodeBoundary[nStart] = std::make_pair(e0, intersection[0]);
				e0 = intersection[0];
				nb_boundary->set(nStart.id(), i);
				nb_boundary->set(nEnd.id(), i);
				nStart = e0.getOppositeNode(nStart);
				listBoundary.push_back(e0);
			}
		}
		auto it = std::find(FrontEdges.begin(), FrontEdges.end(), e0);
		if (it != FrontEdges.end()) {
			// Si l'élément est trouvé, le supprimer
			FrontEdges.erase(it);
		}
		std::get<1>(m_nodeFront.at(nEnd)) = e0;
		std::get<0>(m_nodeBoundary.at(nEnd)) = e0;
		m_listFront.push_back(listBoundary);
		i++;
	}


	std::cout << "============== End Initialisation ==============\n" << std::endl;
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F);
	vtkWriter.setDataOptions(gmds::N | gmds::F);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");
}

/*----------------------------------------------------------------------------*/

Quadfront::STATUS
Quadfront::execute()
{
	return Quadfront::SUCCESS;
}

/*----------------------------------------------------------------------------*/

double
Quadfront::angle(Node Anode)
{
	// Calcul la somme des angles des triangles ayant Anode en commun
	std::vector<Face> vFace = Anode.get<Face>();
	double m_sommeAngle;
	for (auto &face : vFace) {
		if (face.nbEdges() == 3) {
			std::map<Node, std::pair<Edge, Edge>> node2edge = get_edgeInFace(face);
			m_sommeAngle += node2Vector(std::get<0>(node2edge[Anode]), Anode).normalize().angle(node2Vector(std::get<1>(node2edge[Anode]), Anode).normalize());
		}
	}
	return m_sommeAngle;
}

/*----------------------------------------------------------------------------*/

void Quadfront::initialise_Boundary(){
	for (auto& n : m_nodeFront){
		if (angle(n.first)<3*M_PI/4){
			std::get<0>(m_nodeFront[n.first]) = 1;
		}
		else{
			std::get<0>(m_nodeFront[n.first]) = 0;
		}
	}
}

/*----------------------------------------------------------------------------*/
void Quadfront::update_Boundary2(Edge& Aedge_base, Edge& Aedge_up){
	for(auto& node_base : Aedge_base.get<Node>()){
		std::vector<Edge> vEdge;
		bool test_node = false;
		for (auto &e : node_base.get<Edge>()) {
			if (e.get<Face>().size() == 1){
				if (e.get<Face>()[0].nbEdges() == 3){
					vEdge.push_back(e);
					test_node = true;
				}
			}
			else if ((e.get<Face>()[0].nbEdges() == 3 and e.get<Face>()[1].nbEdges() == 4) or
			         e.get<Face>()[0].nbEdges() == 4 and e.get<Face>()[1].nbEdges() == 3) {
				vEdge.push_back(e);
				test_node = true;
			}
		}
		if (test_node) {
			std::get<0>(m_nodeFront[node_base]) = 1;
			std::get<1>(m_nodeFront[node_base]) = vEdge[0];
			std::get<2>(m_nodeFront[node_base]) = vEdge[1];
		}
		else{
			m_nodeFront.erase(m_nodeFront.find(node_base));
		}
	}

	for (auto& node_top : Aedge_up.get<Node>()){
		std::vector<Edge> vEdge;
		bool test_node = false;
		for (auto &e : node_top.get<Edge>()) {
			if (e.get<Face>().size() == 1){
				if (e.get<Face>()[0].nbEdges() == 3){
					vEdge.push_back(e);
					test_node = true;
				}
			}
			else if ((e.get<Face>()[0].nbEdges() == 3 and e.get<Face>()[1].nbEdges() == 4) or
				   (e.get<Face>()[0].nbEdges() == 4 and e.get<Face>()[1].nbEdges() == 3)){
				vEdge.push_back(e);
				test_node = true;
			}
		}
		std::cout<<node_top<<" : "<<vEdge.size()<<std::endl;
		if (test_node) {
			if (m_nodeFront.find(node_top) == m_nodeFront.end()){
				Edge e0;
				m_nodeFront[node_top] = std::make_tuple(0, e0, e0);
			}
			if (vEdge.size() == 2) {
				std::get<1>(m_nodeFront[node_top]) = vEdge[0];
				std::get<2>(m_nodeFront[node_top]) = vEdge[1];
			}
			/*
			else {
				for (auto& ev : vEdge){
					for (auto& la : m_listFront){
						Face quad = ev.get<Face>()[0].nbEdges() == 4 ? ev.get<Face>()[0] : ev.get<Face>()[1];
						Edge edge_opposit = get_oppositEdge_Quad(ev, quad);
						auto it = std::find(la.begin(), la.end(), edge_opposit);
						if (it == la.end()){
							auto x = std::find(vEdge.begin(), vEdge.end(), ev);
							vEdge.erase(x);
						}
					}
				}

				bool test_front = false;
				for (auto& la : m_listFront){
					for (auto& ea : la){
						if (ea.get<Node>()[0] == node_top or ea.get<Node>()[1] == node_top){
							test_front = true;
							break;
						}
					}
				}
				if (test_front){

				}
				Face face_quad = std::get<1>(m_nodeFront[node_top]).get<Face>()[0].nbEdges() == 4 ? std::get<1>(m_nodeFront[node_top]).get<Face>()[0] :
				                                                                                     std::get<1>(m_nodeFront[node_top]).get<Face>()[1];
				Edge edge0 = std::get<1>(m_nodeFront[node_top]);
				Edge edge1 = std::get<2>(m_nodeFront[node_top]);
				for(auto& l : m_listFront){
					for(auto& e_l : l){
						if (e_l == get_oppositEdge_Quad(edge0, face_quad)){
							std::get<2>(m_nodeFront[node_top]) = Aedge_up;
						}
						else if (e_l == get_oppositEdge_Quad(edge1, face_quad)){
							std::get<1>(m_nodeFront[node_top]) = Aedge_up;
						}
					}
				}

			}
			 */
			std::get<0>(m_nodeFront[node_top]) = angle(node_top) < 3*M_PI/4 ? 1 : 0;
		}
		else{
			if (m_nodeFront.find(node_top) != m_nodeFront.end()){
				m_nodeFront.erase(m_nodeFront.find(node_top));
			}
		}
	}
}


/*----------------------------------------------------------------------------*/


void Quadfront::set_nodeBoundary(Node &Anode, int &i, Edge &Aedge0, Edge &Aedge1)
{
	m_nodeFront[Anode] = std::make_tuple(i, Aedge0, Aedge1);
}

/*----------------------------------------------------------------------------*/

std::tuple<int, Edge, Edge> Quadfront::get_nodeBoundary(Node &Anode)
{
	return m_nodeFront[Anode];
}

/*----------------------------------------------------------------------------*/

bool Quadfront::find_nodeBoundary(Node &Anode)
{
	return m_nodeFront.find(Anode) != m_nodeFront.end();
}

/*----------------------------------------------------------------------------*/

std::map<Node, std::pair<Edge, Edge>> Quadfront::get_edgeInFace(Face Aface)
{
	std::map<Node, std::pair<Edge, Edge>> result;

	for (auto& e : Aface.get<Edge>()){
		for(auto& n : e.get<Node>()){
			if (result.find(n) == result.end()){
				result[n] = std::make_pair(e, e);
			}
			else{
				std::get<1>(result[n]) = e;
			}
		}
	}
	return result;
}

/*----------------------------------------------------------------------------*/

Node Quadfront::get_opposit2Edge(Edge &Aedge, Face &Aface)
{
	if (Aface.nbNodes() == 3){
		std::vector<Node> vNode = Aface.get<Node>();
		for (auto &n : vNode) {
			if (n != Aedge.get<Node>()[0] and n != Aedge.get<Node>()[1]) {
				return n;
			}
		}
	}
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::get_opposit2Node(Node &Anode, Face &Aface)
{
	if (Aface.nbNodes() == 3){
		std::vector<Edge> vEdge = Aface.get<Edge>();
		for (auto &e : vEdge) {
			if (Anode != e.template get<Node>()[0] and Anode != e.template get<Node>()[1]) {
				return e;
			}
		}
	}
}

/*----------------------------------------------------------------------------*/

math::Vector3d
Quadfront::node2Vector(Edge &Aedge, Node &Anode)
{
	return {Aedge.getOppositeNode(Anode).X() - Anode.X(), Aedge.getOppositeNode(Anode).Y() - Anode.Y(), Aedge.getOppositeNode(Anode).Z() - Anode.Z()};
}

/*----------------------------------------------------------------------------*/

math::Vector3d Quadfront::get_bissectrice(Node& Anode, Edge& Aedge0, Edge& Aedge1){
	math::Vector3d bissectrice = {0, 0, 0};
	bissectrice = node2Vector(Aedge0, Anode).normalize() + node2Vector(Aedge1, Anode).normalize();

	Node node0 = Aedge0.getOppositeNode(Anode);
	Node node1 = Aedge1.getOppositeNode(Anode);
	math::Vector3d bissectrice1 = node2Vector(Aedge0, node0).normalize() + node2Vector(Aedge1, node1).normalize();
	math::Vector3d vec_norm1 = {bissectrice1[1], -bissectrice1[0], bissectrice1[2]};
	math::Vector3d bissectrice2 = node2Vector(Aedge0, Anode).normalize() + node2Vector(Aedge1, Anode).normalize();
	math::Vector3d vec_norm2 = {bissectrice2[1], -bissectrice2[0], bissectrice2[2]};
	if (bissectrice1.dot(bissectrice) > bissectrice2.dot(bissectrice)){
		return bissectrice1;
	}
	else{
		return bissectrice2;
	}
	//return bissectrice.normalize();
}

/*----------------------------------------------------------------------------*/

math::Vector3d Quadfront::intersectionVec2Edge(math::Vector3d &vec_bissectrice, Node &Anode, Edge &Aedge)
{
	// a1X+b1Y+c1 = 0 pour le vec
	// a2X+b2Y+c2 = 0 pour le edge
	// Anode origine du vecteur vec

	math::Vector3d vecEdge = node2Vector(Aedge, Aedge.get<Node>()[0]);

	double a1 = vec_bissectrice.Y();
	double b1 = -vec_bissectrice.X();
	double c1 = -(a1 * Anode.X() + b1 * Anode.Y());

	double a2 = vecEdge.Y();
	double b2 = -vecEdge.X();
	double c2 = -(a2 * Aedge.get<Node>()[0].X() + b2 * Aedge.get<Node>()[0].Y());

	double y1 = (c1 * a2 / a1) - c2;
	double y2 = (b2 - b1 * a2 / a1);

	double y = y1 / y2;
	double x = -(c1 + b1 * y) / a1;
	return {x, y, 0};
}
/*----------------------------------------------------------------------------*/

void
Quadfront::removeFace(Face &Aface)
{
	for (auto &e : Aface.get<Edge>()) {
		e.remove<Face>(Aface);
		e.get<Node>()[0].remove<Face>(Aface);
		e.get<Node>()[1].remove<Face>(Aface);
	}
	for (auto &n : Aface.get<Node>()) {
		n.remove<Face>(Aface);
	}
	m_mesh->deleteFace(Aface.id());
}

/*----------------------------------------------------------------------------*/

void
Quadfront::removeEdge(Edge &Aedge)
{
	for (auto &f : Aedge.get<Face>()) {
		f.remove<Edge>(Aedge);
	}
	std::vector<Node> vNode = Aedge.get<Node>();
	vNode[0].remove<Edge>(Aedge);
	vNode[1].remove<Edge>(Aedge);
	m_mesh->deleteEdge(Aedge.id());
}

/*----------------------------------------------------------------------------*/

Edge
Quadfront::swap(Node &Anode, Face& Aface, Face& Aface_opposit){

	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(Aface);
	std::map<Node, std::pair<Edge, Edge>> mapFace_opposit = get_edgeInFace(Aface_opposit);

	Edge edge_inter = get_opposit2Node(Anode, Aface);
	Node node_opposit = get_opposit2Edge(edge_inter, Aface_opposit);

	Edge edge0 = std::get<0>(mapFace[Anode]);
	Edge edge1 = std::get<1>(mapFace[Anode]);

	Node node1 = edge0.getOppositeNode(Anode);
	Node node2 = edge1.getOppositeNode(Anode);

	removeFace(Aface);
	removeFace(Aface_opposit);


	// Création du de la diagonale
	Edge sideEdge = m_mesh->newEdge(Anode, node_opposit);

	Anode.add<Edge>(sideEdge);
	node_opposit.add<Edge>(sideEdge);

	// Suppression de l'arete opposée à Anode
	removeEdge(edge_inter);

	// Création des deux nouveaux triangle
	Face newFace1 = m_mesh->newTriangle(Anode, node_opposit, node1);
	Face newFace2 = m_mesh->newTriangle(Anode, node_opposit, node2);

	// mettre à jour les connections des noeuds avec les faces
	Anode.add<Face>(newFace1);
	Anode.add<Face>(newFace2);
	node_opposit.add<Face>(newFace1);
	node_opposit.add<Face>(newFace2);
	node1.add<Face>(newFace1);
	node2.add<Face>(newFace2);

	// mettre à jour la connection des faces avec les edges
	newFace1.add<Edge>(sideEdge);
	newFace1.add<Edge>(edge0);
	newFace2.add<Edge>(sideEdge);
	newFace2.add<Edge>(edge1);

	sideEdge.add<Face>(newFace1);
	edge0.add<Face>(newFace1);
	sideEdge.add<Face>(newFace2);
	edge1.add<Face>(newFace2);

	// On cherche le coté opposé à edge_min
	if (std::get<0>(mapFace_opposit[node1]) != edge_inter) {
		newFace1.add<Edge>(std::get<0>(mapFace_opposit[node1]));
		std::get<0>(mapFace_opposit[node1]).add<Face>(newFace1);
	}
	else {
		newFace1.add<Edge>(std::get<1>(mapFace_opposit[node1]));
		std::get<1>(mapFace_opposit[node1]).add<Face>(newFace1);
	}
	if (std::get<0>(mapFace_opposit[node2]) != edge_inter) {
		newFace2.add<Edge>(std::get<0>(mapFace_opposit[node2]));
		std::get<0>(mapFace_opposit[node2]).add<Face>(newFace2);
	}
	else {
		newFace2.add<Edge>(std::get<1>(mapFace_opposit[node2]));
		std::get<1>(mapFace_opposit[node2]).add<Face>(newFace2);
	}

	return sideEdge;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::split(Node &AnewNode, Node &Anode, Face& Aface){

	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(Aface);
	//std::map<Node, std::pair<Edge, Edge>> mapFace_opposit = get_edgeInFace(Aface_opposit);

	Edge edge_inter = get_opposit2Node(Anode, Aface);
	//Node node_opposit = get_opposit2Edge(edge_inter, Aface_opposit);

	Edge edge1 = std::get<0>(mapFace[Anode]);
	Edge edge2 = std::get<1>(mapFace[Anode]);

	Node node1 = edge1.getOppositeNode(Anode);
	Node node2 = edge2.getOppositeNode(Anode);

	Edge splitEdge1 = AnewNode.get<Edge>()[0].getOppositeNode(AnewNode) == node1 ? AnewNode.get<Edge>()[0] : AnewNode.get<Edge>()[1];
	Edge splitEdge2 = AnewNode.get<Edge>()[1].getOppositeNode(AnewNode) == node2 ? AnewNode.get<Edge>()[1] : AnewNode.get<Edge>()[0];


	Edge sideEdge = m_mesh->newEdge(Anode, AnewNode);

	Anode.add<Edge>(sideEdge);
	AnewNode.add<Edge>(sideEdge);


	// Création des deux nouveaux triangle
	Face newFace1 = m_mesh->newTriangle(Anode, AnewNode, node1);
	Face newFace2 = m_mesh->newTriangle(Anode, AnewNode, node2);

	// mettre à jour les connections des noeuds avec les faces
	Anode.add<Face>(newFace1);
	Anode.add<Face>(newFace2);
	AnewNode.add<Face>(newFace1);
	AnewNode.add<Face>(newFace2);
	node1.add<Face>(newFace1);
	node2.add<Face>(newFace2);


	// mettre à jour la connection des faces avec les edges
	// CHECKER
	newFace1.add<Edge>(sideEdge);
	newFace2.add<Edge>(sideEdge);

	newFace1.add<Edge>(splitEdge1);
	newFace2.add<Edge>(splitEdge2);

	newFace1.add<Edge>(edge1);
	newFace2.add<Edge>(edge2);

	sideEdge.add<Face>(newFace1);
	sideEdge.add<Face>(newFace2);

	splitEdge1.add<Face>(newFace1);
	splitEdge2.add<Face>(newFace2);
	edge1.add<Face>(newFace1);
	edge2.add<Face>(newFace2);


	// Supprimer les faces
	removeFace(Aface);

	return sideEdge;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::split_operation(Node &AnewNode, Node &Anode, Face& Aface, Face& Aface_opposit){
	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(Aface);

	Edge edge_inter = get_opposit2Node(Anode, Aface);
	Node node_opposit = get_opposit2Edge(edge_inter, Aface_opposit);

	Edge edge1 = std::get<0>(mapFace[Anode]);
	Edge edge2 = std::get<1>(mapFace[Anode]);

	Node node1 = edge1.getOppositeNode(Anode);
	Node node2 = edge2.getOppositeNode(Anode);

	Edge splitEdge1 = m_mesh->newEdge(node1, AnewNode);
	Edge splitEdge2 = m_mesh->newEdge(AnewNode, node2);

	AnewNode.add<Edge>(splitEdge1);
	node1.add<Edge>(splitEdge1);
	AnewNode.add<Edge>(splitEdge2);
	node2.add<Edge>(splitEdge2);

	Edge sideEdge1 = split(AnewNode, Anode, Aface);
	Edge sideEdge2 = split(AnewNode, node_opposit, Aface_opposit);

	removeEdge(edge_inter);

	return sideEdge1;
}

/*----------------------------------------------------------------------------*/

// checker au préalable le statut du noeud => dans le cas ou c'est 0
// < pi/6 pour que ce soit valable
Edge Quadfront::sideEdge(Node &Anode, Edge &Aedge)
{

	Edge sideEdge;


	// Bissectrice des arêtes adjacentes au noeud
	math::Vector3d m_bissectrice = {0, 0, 0};
	for (auto& f : Anode.get<Face>()){
		if (f.nbNodes() == 3) {
			std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(f);
			m_bissectrice += get_bissectrice(Anode, std::get<0>(mapFace[Anode]), std::get<1>(mapFace[Anode])).normalize();
		}
	}
	m_bissectrice.normalize();

	// Obtention du meilleur edge

	double dot_max1;
	double dot_max2;
	Edge edge_min1;
	Edge edge_min2;

	std::vector<Edge> vEdge = Anode.get<Edge>();
	std::vector<Edge> vEdge_bis;


	for (auto ed : vEdge){
		if (ed.get<Face>().size() == 1){
			vEdge_bis.push_back(ed);
		}
		else if (ed.get<Face>()[0].nbEdges() == 3 or ed.get<Face>()[1].nbEdges() == 3){
			vEdge_bis.push_back(ed);
		}
	}

	vEdge = vEdge_bis;
	std::cout<<vEdge.size()<<std::endl;

	if (node2Vector(vEdge[0], Anode).normalize().dot(m_bissectrice) > node2Vector(vEdge[1], Anode).normalize().dot(m_bissectrice)) {
		dot_max1 = node2Vector(vEdge[0], Anode).normalize().dot(m_bissectrice);
		dot_max2 = node2Vector(vEdge[1], Anode).normalize().dot(m_bissectrice);
		edge_min1 = vEdge[0];
		edge_min2 = vEdge[1];
	}
	else {
		dot_max1 = node2Vector(vEdge[1], Anode).normalize().dot(m_bissectrice);
		dot_max2 = node2Vector(vEdge[0], Anode).normalize().dot(m_bissectrice);
		edge_min1 = vEdge[1];
		edge_min2 = vEdge[0];
	}

	std::cout << "Edge " << edge_min1.id() << " : " << edge_min1.get<Node>()[0] << " ; " << edge_min1.get<Node>()[1] << " : " << dot_max1 << std::endl;
	std::cout << "Edge " << edge_min2.id() << " : " << edge_min2.get<Node>()[0] << " ; " << edge_min2.get<Node>()[1] << " : " << dot_max2 << std::endl;

	for (Edge &e : vEdge) {
		if (e != edge_min1 and e != edge_min2) {
			math::Vector3d vec = node2Vector(e, Anode).normalize();
			std::cout << "Edge " << e.id() << " : " << e.get<Node>()[0] << " ; " << e.get<Node>()[1] << " : " << m_bissectrice.dot(vec) << std::endl;
			if (m_bissectrice.dot(vec) > dot_max1){
				dot_max2 = dot_max1;
				dot_max1 = m_bissectrice.dot(vec);
				edge_min2 = edge_min1;
				edge_min1 = e;
			}
			else if (m_bissectrice.dot(vec) > dot_max2) {
				dot_max2 = m_bissectrice.dot(vec);
				edge_min2 = e;
			}
		}
	}

	std::cout << "edge_min1 " << edge_min1.id() << " : " << edge_min1.get<Node>()[0] << " ; " << edge_min1.get<Node>()[1] << " : " << dot_max1 << std::endl;
	std::cout << "edge_min2 " << edge_min2.id() << " : " << edge_min2.get<Node>()[0] << " ; " << edge_min2.get<Node>()[1] << " : " << dot_max2 << std::endl;

	Edge Aedge2 = std::get<1>(m_nodeFront[Anode]) != Aedge ? std::get<1>(m_nodeFront[Anode]) : std::get<2>(m_nodeFront[Anode]);
	Edge edge0 = std::get<1>(m_nodeFront[Aedge.getOppositeNode(Anode)]) != Aedge ? std::get<1>(m_nodeFront[Aedge.getOppositeNode(Anode)]) :
	                                                                                std::get<2>(m_nodeFront[Aedge.getOppositeNode(Anode)]);
	Edge edge1 = std::get<1>(m_nodeFront[Aedge2.getOppositeNode(Anode)]) != Aedge2 ? std::get<1>(m_nodeFront[Aedge2.getOppositeNode(Anode)]) :
	                                                                                std::get<2>(m_nodeFront[Aedge2.getOppositeNode(Anode)]);


	// CAS IDEAL

	if (acos(dot_max1) < M_PI / 6) {
		Node node_opposit = edge_min1.getOppositeNode(Anode);
		if (m_nodeFront.find(node_opposit) == m_nodeFront.end()){
			sideEdge = edge_min1;
			return sideEdge;
		}
		if (node_opposit != edge0.getOppositeNode(Aedge.getOppositeNode(Anode)) and node_opposit != edge1.getOppositeNode(Aedge2.getOppositeNode(Anode))){
			sideEdge = edge_min1;
			std::cout << "l'angle " << edge_min1.id() << " est plus petit que Pi/12 " << std::endl;
			return sideEdge;
		}
	}
	//bool testb = edge_min1.getOppositeNode(Anode) != edge0.getOppositeNode(Aedge.getOppositeNode(Anode));
	std::cout<<"testb : "<<edge0.id()<<std::endl;

	if (edge_min1.getOppositeNode(Anode) != edge0.getOppositeNode(Aedge.getOppositeNode(Anode)) and
	    edge_min1.getOppositeNode(Anode) != edge1.getOppositeNode(Aedge2.getOppositeNode(Anode)) and
	    m_nodeFront.find(edge_min1.getOppositeNode(Anode)) != m_nodeFront.end() and edge_min1 != Aedge){
		sideEdge = edge_min1;
		std::cout << "l'angle " << edge_min1.id() << " est plus petit que Pi/12 " << std::endl;
		return sideEdge;
	}

	std::cout<<"ici"<<std::endl;

	if (acos(dot_max2) < M_PI / 6) {
		Node node_opposit = edge_min2.getOppositeNode(Anode);
		if (m_nodeFront.find(node_opposit) == m_nodeFront.end()){
			sideEdge = edge_min2;
			return sideEdge;
		}
		if (node_opposit != edge0.getOppositeNode(Aedge.getOppositeNode(Anode)) and node_opposit != edge1.getOppositeNode(Aedge2.getOppositeNode(Anode))) {
			sideEdge = edge_min2;
			return sideEdge;
		}
	}


	if (edge_min2.getOppositeNode(Anode) != edge0.getOppositeNode(Aedge.getOppositeNode(Anode)) and
	    edge_min2.getOppositeNode(Anode) != edge1.getOppositeNode(Aedge2.getOppositeNode(Anode)) and
	    m_nodeFront.find(edge_min2.getOppositeNode(Anode)) != m_nodeFront.end() and edge_min2 != Aedge){
		sideEdge = edge_min2;
		return sideEdge;
	}


	if (vEdge.size() == 3) {

		if (Aedge.getOppositeNode(Anode).get<Edge>().size() == 4) {
			edge_min2 = std::get<1>(m_nodeFront[Anode]) == Aedge ? std::get<2>(m_nodeFront[Anode]) : std::get<1>(m_nodeFront[Anode]);
		}
		else {
			math::Vector3d vec_edge0 = node2Vector(Aedge, Anode).normalize();
			math::Vector3d vecNormal_edge0 = {vec_edge0[1], -vec_edge0[0], vec_edge0[2]};
			math::Vector3d vec_edgemin = node2Vector(edge_min1, Anode).normalize();
			math::Vector3d vecNormal_edgemin = {vec_edgemin[1], -vec_edgemin[0], vec_edgemin[2]};
			if ((vecNormal_edgemin.dot(m_bissectrice) > 0 and vecNormal_edge0.dot(m_bissectrice) < 0)
				 or (vecNormal_edgemin.dot(m_bissectrice) < 0 and vecNormal_edge0.dot(m_bissectrice) > 0)) {
				edge_min2 = Aedge;
			}
			else {
				edge_min2 = std::get<2>(m_nodeFront[Anode]);
			}
		}
	}

	math::Vector3d vec_edgemin1 = node2Vector(edge_min1, Anode).normalize();
	math::Vector3d vecNormal_edgemin1 = {vec_edgemin1[1], -vec_edgemin1[0], vec_edgemin1[2]};
	math::Vector3d vec_edgemin2 = node2Vector(edge_min2, Anode).normalize();
	math::Vector3d vecNormal_edgemin2 = {vec_edgemin2[1], -vec_edgemin2[0], vec_edgemin2[2]};
	if ((vecNormal_edgemin1.dot(m_bissectrice) > 0 and vecNormal_edgemin2.dot(m_bissectrice) > 0)
		 or (vecNormal_edgemin1.dot(m_bissectrice) < 0 and vecNormal_edgemin2.dot(m_bissectrice) < 0)) {
		m_bissectrice = (vec_edgemin1 + vec_edgemin2).normalize() ;
	}

	Node node1 = edge_min1.getOppositeNode(Anode);
	Node node2 = edge_min2.getOppositeNode(Anode);

	std::vector<TCellID> vFace = m_mesh->getCommonFaces(node1, node2);

	bool test = false;
	Face face;
	Face face_opposit;

	for (auto& n : m_mesh->get<Face>(vFace[0]).get<Node>()){
		if (n == Anode){
			test = true;
			face = m_mesh->get<Face>(vFace[0]);
			face_opposit = m_mesh->get<Face>(vFace[1]);
		}
	}
	if (!test){
		face = m_mesh->get<Face>(vFace[1]);
		face_opposit = m_mesh->get<Face>(vFace[0]);
	}

	Edge edge_inter = get_opposit2Node(Anode, face);
	if (std::get<1>(m_nodeFront[Aedge2.getOppositeNode(Anode)]) == edge_inter or std::get<2>(m_nodeFront[Aedge2.getOppositeNode(Anode)]) == edge_inter){
		return Aedge2;
	}
	Node node_opposit = get_opposit2Edge(edge_inter, face_opposit);

	math::Vector3d vec = {node_opposit.X() - Anode.X(), node_opposit.Y() - Anode.Y(), node_opposit.Z() - Anode.Z()};


	// SWAP CASE
	if (m_bissectrice.angle(vec) < M_PI / 6 and vec.norm() <  sqrt(3)*(edge_min1.length() + edge_min2.length()) / 4 ) {
		std::cout << "==== SWAP ====" << std::endl;
		sideEdge = swap(Anode, face, face_opposit);
	}

	// SPLIT CASE
	else {
		std::cout << "edge_inter : "<<edge_inter.get<Node>()[0]<<" ; "<<edge_inter.get<Node>()[1]<<std::endl;
		if (m_nodeFront.find(edge_inter.get<Node>()[0]) != m_nodeFront.end() and m_nodeFront.find(edge_inter.get<Node>()[1]) != m_nodeFront.end()){
			if (std::get<1>(m_nodeFront[edge_inter.get<Node>()[0]]) == edge_inter or std::get<2>(m_nodeFront[edge_inter.get<Node>()[1]]) == edge_inter){

			}
		}
		math::Vector3d vNode = intersectionVec2Edge(m_bissectrice, Anode, edge_inter);
		Node newNode = m_mesh->newNode(vNode.X(), vNode.Y(), vNode.Z());

		std::cout << "==== SPLIT ====" << std::endl;
		Node node_opposit = get_opposit2Edge(edge_inter, face_opposit);
		sideEdge = split_operation(newNode, Anode, face, face_opposit);
	}

	return sideEdge;
}
/*----------------------------------------------------------------------------*/

std::vector<Edge> Quadfront::interestionS(Node &AnodeRecovery0, Node &AnodeRecovery1){

	math::Vector3d vecAB({AnodeRecovery1.X() - AnodeRecovery0.X(), AnodeRecovery1.Y() - AnodeRecovery0.Y(), AnodeRecovery1.Z() - AnodeRecovery0.Z()});

	math::Vector3d vecNormaleAB({vecAB[1], -vecAB[0], vecAB[2]});

	// CHECKER POUR LA 3D
	std::vector<Face> vFace0 = AnodeRecovery0.get<Face>();
	std::vector<Edge> vEdgeIntersect;
	std::map<Node, std::pair<Edge, Edge>> node2edgeInFace;
	Edge e_opposit;
	Face face;

	for (auto& e : AnodeRecovery0.get<Edge>()){
		if (e.getOppositeNode(AnodeRecovery0) == AnodeRecovery1){
			vEdgeIntersect.push_back(e);
			return vEdgeIntersect;
		}
	}


	// determine la meilleure face (celle dont une arete se situe en dessous et l'autre au dessus)
	for (auto &f : vFace0) {

		node2edgeInFace = get_edgeInFace(f);

		math::Vector3d vec0 = node2Vector(std::get<0>(node2edgeInFace[AnodeRecovery0]), AnodeRecovery0);
		math::Vector3d vecNormal0 = {vec0[1], -vec0[0], vec0[2]};
		math::Vector3d vec1 = node2Vector(std::get<1>(node2edgeInFace[AnodeRecovery0]), AnodeRecovery0);
		math::Vector3d vecNormal1 = {vec1[1], -vec1[0], vec0[2]};

		// Dans le cas ou une arête est alignée avec le segment AB
		if (vecAB.dot(vec0) == 1) {
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).X() += 0.05 * vecNormaleAB.X();
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Y() += 0.05 * vecNormaleAB.Y();
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Z() += 0.05 * vecNormaleAB.Z();
			vec0 = vec0 + 0.05 * vecNormaleAB;
			vecNormal0 = {vec0[1], -vec0[0], vec0[2]};
		}
		if (vecAB.dot(vec1) == 1) {
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).X() += 0.05 * vecNormaleAB.X();
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Y() += 0.05 * vecNormaleAB.Y();
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Z() += 0.05 * vecNormaleAB.Z();
			vec1 = vec1 + 0.05 * vecNormaleAB;
			vecNormal1 = {vec1[1], -vec1[0], vec1[2]};
		}

		if ((vecAB.dot(vecNormal0) > 0 and vecAB.dot(vecNormal1) < 0) or (vecAB.dot(vecNormal0) < 0 and vecAB.dot(vecNormal1) > 0)) {
			math::Vector3d vec_bis = (vec0.normalize() + vec1.normalize()).normalize();
			if (vec_bis.dot(vecAB) > 0) {
				e_opposit = get_opposit2Node(AnodeRecovery0, f);
				face = f;
			}
		}
	}

	int i = 0;

	std::cout<<"e_opposit : "<<e_opposit.get<Node>()[0]<<" ; "<<e_opposit.get<Node>()[1]<<std::endl;

	bool test = true;

	if (m_nodeFront.find(e_opposit.get<Node>()[0]) != m_nodeFront.end()){
		if (std::get<1>(m_nodeFront[e_opposit.get<Node>()[0]]) == e_opposit or std::get<2>(m_nodeFront[e_opposit.get<Node>()[0]]) == e_opposit){
			test = false;
		}
	}
	if (m_nodeFront.find(e_opposit.get<Node>()[1]) != m_nodeFront.end()){
		if (std::get<1>(m_nodeFront[e_opposit.get<Node>()[1]]) == e_opposit or std::get<2>(m_nodeFront[e_opposit.get<Node>()[1]]) == e_opposit){
			test = false;
		}
	}

	while (test) {
		i++;


		vEdgeIntersect.push_back(e_opposit);

		face = (e_opposit.get<Face>()[0] == face) ? e_opposit.get<Face>()[1] : e_opposit.get<Face>()[0];

		node2edgeInFace = get_edgeInFace(face);

		if (node2edgeInFace.find(AnodeRecovery1) == node2edgeInFace.end()) {
			Node n_opposit = get_opposit2Edge(e_opposit, face);

			math::Vector3d vecAC = {n_opposit.X() - AnodeRecovery0.X(), n_opposit.Y() - AnodeRecovery0.Y(), n_opposit.Z() - AnodeRecovery0.Z()};

			double test = 1 - vecAB.dot(vecAC);

			if (test == 0.) {
				if (m_nodeFront.find(n_opposit) == m_nodeFront.end()) {
					n_opposit.X() += 0.05 * vecNormaleAB.X();
					n_opposit.Y() += 0.05 * vecNormaleAB.Y();
					n_opposit.Z() += 0.05 * vecNormaleAB.Z();
					vecAC = vecAC + 0.05 * vecNormaleAB;
				}
			}
			math::Vector3d vecNormaleAC = {vecAC[1], -vecAC[0], vecAC[2]};

			Node node0 = std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit);

			math::Vector3d vec0 = {std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).X() - AnodeRecovery0.X(),
			                       std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Y() - AnodeRecovery0.Y(),
			                       std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Z() - AnodeRecovery0.Z()};

			math::Vector3d vecNormale0 = {vec0[1], -vec0[0], vec0[2]};

			math::Vector3d vec1 = {std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).X() - AnodeRecovery0.X(),
			                       std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Y() - AnodeRecovery0.Y(),
			                       std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Z() - AnodeRecovery0.Z()};

			math::Vector3d vecNormale1 = {vec0[1], -vec0[0], vec0[2]};

			if (vecAB.dot(vecNormaleAC) > 0) {
				e_opposit = (vecAB.dot(vecNormale0) > 0) ? std::get<1>(node2edgeInFace[n_opposit]) : std::get<0>(node2edgeInFace[n_opposit]);
			}
			else {
				e_opposit = (vecAB.dot(vecNormale0) > 0) ? std::get<0>(node2edgeInFace[n_opposit]) : e_opposit = std::get<1>(node2edgeInFace[n_opposit]);
			}
		}
		else {
			break;
		}

		if (m_nodeFront.find(e_opposit.get<Node>()[0]) != m_nodeFront.end()){
			if (std::get<1>(m_nodeFront[e_opposit.get<Node>()[0]]) == e_opposit or std::get<2>(m_nodeFront[e_opposit.get<Node>()[0]]) == e_opposit){
				test = false;
			}
		}
		if (m_nodeFront.find(e_opposit.get<Node>()[1]) != m_nodeFront.end()){
			if (std::get<1>(m_nodeFront[e_opposit.get<Node>()[1]]) == e_opposit or std::get<2>(m_nodeFront[e_opposit.get<Node>()[1]]) == e_opposit){
				test = false;
			}
		}
	}
	return vEdgeIntersect;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::recoveryEdge(std::vector<Edge>  listInside, Node& nodeRecovery0, Node& nodeRecovery1){

	if (listInside.size() == 1){
		if (listInside[0].get<Node>()[0] == nodeRecovery0 or listInside[0].get<Node>()[0] == nodeRecovery1){
			return listInside[0];
		}
	}

	Edge edgeRecovery;

	math::Vector3d vecAB({nodeRecovery1.X()-nodeRecovery0.X(),
	                      nodeRecovery1.Y()-nodeRecovery0.Y(),
	                      nodeRecovery1.Z()-nodeRecovery0.Z()});

	for (auto e = listInside.begin(); e!=listInside.end();){
		std::vector<Node> vNode = (*e).get<Node>();
		std::vector<Face> vFace = (*e).get<Face>();
		Node n_opposit0 = get_opposit2Edge(*e, vFace[0]);
		std::map<Node, std::pair<Edge, Edge>> node2edgeInFace0 = get_edgeInFace(vFace[0]);
		Node n_opposit1 = get_opposit2Edge(*e, vFace[1]);
		std::map<Node, std::pair<Edge, Edge>> node2edgeInFace1 = get_edgeInFace(vFace[1]);
		math::Triangle triangle0(n_opposit0.point(), n_opposit1.point(), vNode[1].point());
		math::Triangle triangle1(n_opposit0.point(), n_opposit1.point(), vNode[0].point());


		if (triangle0.area()>0 and triangle1.area()>0){
			Edge diagEdge0;
			bool test_e = false;

			for (auto& k : (*e).get<Face>()[0].get<Node>()){
				if (k == n_opposit0){
					test_e = true;
					diagEdge0 = swap(n_opposit0,  (*e).get<Face>()[0], (*e).get<Face>()[1]);
					break;
				}
			}

			if (!test_e){
				diagEdge0 = swap(n_opposit0,  (*e).get<Face>()[1], (*e).get<Face>()[0]);
			}


			if ((diagEdge0.get<Node>()[0] == nodeRecovery0 and diagEdge0.get<Node>()[1] == nodeRecovery1) or
			    (diagEdge0.get<Node>()[1] == nodeRecovery0 and diagEdge0.get<Node>()[0] == nodeRecovery1)){
				edgeRecovery = diagEdge0;
			}

			listInside.erase(e);

			math::Vector3d vecNode = intersectionVec2Edge(vecAB, nodeRecovery0, diagEdge0);
			math::Point point = (vecNode.X(), vecNode.Y(), vecNode.Z());
			if (diagEdge0.segment().isIn(point)){
				listInside.push_back(diagEdge0);
			}
		}
		else{
			listInside.erase(e);
			listInside.push_back(*e);
		}
	}
	return edgeRecovery;
}

/*----------------------------------------------------------------------------*/

void Quadfront::testInside(Face &faceRef, Face &faceInside, std::vector<Face> &listInside){
	std::vector<Edge> edgeFaceRef = faceRef.get<Edge>();

	if (std::find(listInside.begin(), listInside.end(), faceInside) != listInside.end() or faceInside.nbEdges() == 4) {
		return;
	}
	bool test = true;
	double Xmax = std::max(std::max(faceRef.get<Node>()[0].X(), faceRef.get<Node>()[1].X()), std::max(faceRef.get<Node>()[2].X(), faceRef.get<Node>()[3].X()));
	double Xmin = std::min(std::min(faceRef.get<Node>()[0].X(), faceRef.get<Node>()[1].X()), std::min(faceRef.get<Node>()[2].X(), faceRef.get<Node>()[3].X()));
	double Ymax = std::max(std::max(faceRef.get<Node>()[0].Y(), faceRef.get<Node>()[1].Y()), std::max(faceRef.get<Node>()[2].Y(), faceRef.get<Node>()[3].Y()));
	double Ymin = std::min(std::min(faceRef.get<Node>()[0].Y(), faceRef.get<Node>()[1].Y()), std::min(faceRef.get<Node>()[2].Y(), faceRef.get<Node>()[3].Y()));
	double Zmax = std::max(std::max(faceRef.get<Node>()[0].Z(), faceRef.get<Node>()[1].Z()), std::max(faceRef.get<Node>()[2].Z(), faceRef.get<Node>()[3].Z()));
	double Zmin = std::min(std::min(faceRef.get<Node>()[0].Z(), faceRef.get<Node>()[1].Z()), std::min(faceRef.get<Node>()[2].Z(), faceRef.get<Node>()[3].Z()));
	for (auto &n : faceInside.get<Node>()) {
		if ((n.X() < Xmin or n.X() > Xmax) or (n.Y() < Ymin or n.Y() > Ymax) or (n.Z() < Zmin or n.Z() > Zmax)) {
			test = false;
		}
	}
	if (test) {
		listInside.push_back(faceInside);
		for (auto &e : faceInside.get<Edge>()) {
			if (std::find(edgeFaceRef.begin(), edgeFaceRef.end(), e) == edgeFaceRef.end()) {
				if (e.get<Face>().size() == 2) {
					if (e.get<Face>()[0] != faceInside) {
						testInside(faceRef, e.get<Face>()[0], listInside);
					}
					else {
						testInside(faceRef, e.get<Face>()[1], listInside);
					}
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/

Face Quadfront::createQuad(Node &Anode, Edge &Aedge, Edge &sideEdge0, Edge &sideEdge1, Edge &recoveryEdge){
	Face face0 = (Aedge.get<Face>()[0].nbEdges() == 3) ? Aedge.get<Face>()[0] : Aedge.get<Face>()[1];

	std::vector<Face> listInside;
	Face quad;

	Node node0 = Aedge.getOppositeNode(Anode);

	if (sideEdge0.get<Node>()[0] == Anode or sideEdge0.get<Node>()[1] == Anode) {
		quad = m_mesh->newQuad(Anode, sideEdge0.getOppositeNode(Anode), sideEdge1.getOppositeNode(node0), node0);
		Anode.add<Face>(quad);
		sideEdge0.getOppositeNode(Anode).add<Face>(quad);
		sideEdge1.getOppositeNode(node0).add<Face>(quad);
		node0.add<Face>(quad);
	}
	else {
		quad = m_mesh->newQuad(Anode, sideEdge1.getOppositeNode(Anode), sideEdge0.getOppositeNode(node0), node0);
		Anode.add<Face>(quad);
		sideEdge1.getOppositeNode(Anode).add<Face>(quad);
		sideEdge0.getOppositeNode(node0).add<Face>(quad);
		node0.add<Face>(quad);
	}


	quad.add<Edge>(Aedge);
	quad.add<Edge>(sideEdge0);
	quad.add<Edge>(sideEdge1);
	quad.add<Edge>(recoveryEdge);

	testInside(quad, face0, listInside);

	Aedge.add<Face>(quad);
	sideEdge0.add<Face>(quad);
	sideEdge1.add<Face>(quad);
	recoveryEdge.add<Face>(quad);


	for (auto &f : listInside) {
		for (auto &e : f.get<Edge>()) {
			if (e != Aedge and e != sideEdge0 and e != sideEdge1 and e != recoveryEdge) {
				removeEdge(e);
			}
		}
		removeFace(f);
	}

	return quad;
}

/*----------------------------------------------------------------------------*/

Node Quadfront::get_oppositNode_Quad(Node &Anode, Face &Aface){
	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(Aface);
	Edge edge = std::get<0>(mapFace[Anode]);
	if (std::get<0>(mapFace[edge.getOppositeNode(Anode)]) == edge){
		return std::get<1>(mapFace[edge.getOppositeNode(Anode)]).getOppositeNode(edge.getOppositeNode(Anode));
	}
	else{
		return std::get<0>(mapFace[edge.getOppositeNode(Anode)]).getOppositeNode(edge.getOppositeNode(Anode));
	}
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::get_oppositEdge_Quad(Edge &Aedge, Face &Aface){
	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(Aface);
	Node node0 = Aedge.get<Node>()[0];
	Edge edge = std::get<0>(mapFace[node0]) == Aedge ? std::get<1>(mapFace[node0]) : std::get<0>(mapFace[node0]);
	if (std::get<0>(mapFace[edge.getOppositeNode(node0)]) == edge){
		return std::get<1>(mapFace[edge.getOppositeNode(node0)]);
	}
	else{
		return std::get<0>(mapFace[edge.getOppositeNode(node0)]);
	}
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::closingSeam(Edge& Aedge0, Edge& Aedge1, Edge& Aedge_triangle0, Edge& Aedge_triangle1){
	Node node_common = Aedge0.get<Node>()[0] == Aedge1.get<Node>()[0] or
	   Aedge0.get<Node>()[0] == Aedge1.get<Node>()[1] ? Aedge0.get<Node>()[0] : Aedge0.get<Node>()[1];

	Node node0 = Aedge0.getOppositeNode(node_common);
	Node node1 = Aedge1.getOppositeNode(node_common);

	Face triangle0 = Aedge_triangle0.get<Face>()[0];
	Face triangle1 = Aedge_triangle1.get<Face>()[0];

	Face quad0 = Aedge0.get<Face>()[0].nbEdges() == 4 ? Aedge0.get<Face>()[0] : Aedge0.get<Face>()[1];
	Face quad1 = Aedge1.get<Face>()[0].nbEdges() == 4 ? Aedge1.get<Face>()[0] : Aedge1.get<Face>()[1];

	std::map<Node, std::pair<Edge, Edge>> mapFace_quad1 = get_edgeInFace(quad1);
	Edge side_delete = std::get<0>(mapFace_quad1[node1]) == Aedge1 ? std::get<1>(mapFace_quad1[node1]) : std::get<0>(mapFace_quad1[node1]);

	Face quad2 = side_delete.get<Face>()[0] == quad1 ? side_delete.get<Face>()[1] : side_delete.get<Face>()[0];
	Edge next_Aedge1 = std::get<1>(m_nodeFront[node1]) == Aedge1 ? std::get<2>(m_nodeFront[node1]) : std::get<1>(m_nodeFront[node1]);

	node_common.remove<Edge>(Aedge1);

	quad1.remove<Edge>(side_delete);
	quad1.remove<Edge>(Aedge1);
	quad1.remove<Node>(node1);
	quad2.remove<Edge>(next_Aedge1);
	quad2.remove<Node>(node1);

	quad1.add<Node>(node0);
	quad2.add<Node>(node0);

	Edge side_quad;
	if (next_Aedge1 != side_delete) {
		quad2.remove<Edge>(side_delete);
		side_quad = m_mesh->newEdge(node0, side_delete.getOppositeNode(node1));
		quad2.add<Edge>(side_quad);
		quad1.add<Edge>(side_quad);
		side_quad.add<Face>(quad1);
		side_quad.add<Face>(quad2);
	}

	Edge new_top2 = m_mesh->newEdge(node0, next_Aedge1.getOppositeNode(node1));

	quad1.add<Edge>(Aedge0);
	//quad1.add<Edge>(side_quad);
	quad2.add<Edge>(new_top2);

	node0.add<Face>(quad1);
	node0.add<Face>(quad2);

	new_top2.add<Face>(quad2);
	Aedge0.add<Face>(quad1);

	triangle1.remove<Node>(node1);
	triangle1.remove<Edge>(Aedge_triangle1);
	triangle1.remove<Edge>(next_Aedge1);

	triangle1.add<Node>(node0);
	triangle1.add<Edge>(new_top2);
	triangle1.add<Edge>(Aedge_triangle0);

	new_top2.add<Face>(triangle1);
	Aedge_triangle0.add<Face>(triangle1);

	removeEdge(Aedge1);
	removeEdge(side_delete);
	removeEdge(next_Aedge1);
	removeEdge(Aedge_triangle1);

	for (auto& e_node1 : node1.get<Edge>()){
		removeEdge(e_node1);
	}

	Edge quad0_base = get_oppositEdge_Quad(Aedge0, quad0);
	update_Boundary2(Aedge0, quad0_base);

	return new_top2;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::transitionSeam(Node& Anode, Edge& Aedge_max, Edge& Aedge_min){

	math::Point center_seg = Aedge_max.center();
	Node node_new = m_mesh->newNode(center_seg.X(), center_seg.Y(), center_seg.Z());
	Edge Aedge_new0 = m_mesh->newEdge(node_new, Anode);
	Edge Aedge_new1 = m_mesh->newEdge(Aedge_max.getOppositeNode(Anode), node_new);

	node_new.add<Edge>(Aedge_new0);
	Anode.add<Edge>(Aedge_new0);
	node_new.add<Edge>(Aedge_new1);
	Aedge_max.getOppositeNode(Anode).add<Edge>(Aedge_new1);

	Face triangle_big = Aedge_max.get<Face>()[0].nbEdges() == 4 ? Aedge_max.get<Face>()[1] : Aedge_max.get<Face>()[0];
	Node nodeTop = get_opposit2Edge(Aedge_max, triangle_big);

	Edge edge_split = split(node_new, nodeTop, triangle_big);


	Face quadBig = Aedge_max.get<Face>()[0].nbEdges() == 4 ? Aedge_max.get<Face>()[0] : Aedge_max.get<Face>()[1];
	std::map<Node, std::pair<Edge, Edge>> mapFace_quadBig = get_edgeInFace(quadBig);
	Edge edge_side = std::get<0>(mapFace_quadBig[Anode]) == Aedge_max ?
	                                                                  std::get<1>(mapFace_quadBig[Anode]) :
	                                                                  std::get<0>(mapFace_quadBig[Anode]);

	Node node_opposit = edge_side.getOppositeNode(Anode);
	Node node_top = Aedge_min.getOppositeNode(Anode);

	Edge edge_new = m_mesh->newEdge(node_new, node_opposit);
	node_new.add<Edge>(edge_new);
	node_opposit.add<Edge>(edge_new);

	quadBig.remove<Node>(Anode);
	Anode.remove<Face>(quadBig);
	quadBig.remove<Edge>(edge_side);
	edge_side.remove<Face>(quadBig);
	quadBig.remove<Edge>(Aedge_max);
	Aedge_max.remove<Face>(quadBig);

	quadBig.add<Node>(node_new);
	node_new.add<Face>(quadBig);
	quadBig.add<Edge>(edge_new);
	edge_new.add<Face>(quadBig);
	quadBig.add<Edge>(Aedge_new1);
	Aedge_new1.add<Face>(quadBig);

	std::vector<Face> listInside;
	std::vector<Edge> listIntersect = interestionS(node_new, node_top);
	Edge edgeRecovery = recoveryEdge(listIntersect, node_new, node_top);

	Face quad_new = createQuad(node_new, edgeRecovery, edge_new, Aedge_min, edge_side);

	removeEdge(Aedge_max);

	return edgeRecovery;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::transitionSplit(Node& Anode, Edge& Aedge_max, Edge& Aedge_min){
	math::Point center_seg = Aedge_max.center();
	Node node_max = Aedge_max.getOppositeNode(Anode);
	Node node_min = Aedge_min.getOppositeNode(Anode);

	Face quadBig = Aedge_max.get<Face>()[0].nbEdges() == 4 ? Aedge_max.get<Face>()[0] : Aedge_max.get<Face>()[1];
	std::map<Node, std::pair<Edge, Edge>> mapFace_quadBig = get_edgeInFace(quadBig);
	std::cout<<"Anode "<<Anode<<std::endl;
	std::cout<<"Node max "<<node_max<<std::endl;
	Edge edge_side;
	if (std::get<0>(mapFace_quadBig[node_max]) == Aedge_max){
		edge_side = std::get<1>(mapFace_quadBig[node_max]);
	}
	else{
		edge_side = std::get<0>(mapFace_quadBig[node_max]);
	}

	std::cout<<"edge_side : "<<edge_side.get<Node>()[0]<<" ; "<<edge_side.get<Node>()[1]<<std::endl;


	Node node_new = m_mesh->newNode(center_seg.X(), center_seg.Y(), center_seg.Z());
	Edge Aedge_new0 = m_mesh->newEdge(node_new, Anode);
	Edge Aedge_new1 = m_mesh->newEdge(node_max, node_new);

	node_new.add<Edge>(Aedge_new0);
	Anode.add<Edge>(Aedge_new0);
	node_new.add<Edge>(Aedge_new1);
	node_max.add<Edge>(Aedge_new1);

	Face triangle_big = Aedge_max.get<Face>()[0].nbEdges() == 4 ? Aedge_max.get<Face>()[1] : Aedge_max.get<Face>()[0];
	Node nodeTop = get_opposit2Edge(Aedge_max, triangle_big);

	Edge edge_split = split(node_new, nodeTop, triangle_big);

	Node node_opposit = get_oppositNode_Quad(Anode, quadBig);

	math::Point center_diag = quadBig.center();
	Node nodeQuad_split = m_mesh->newNode(center_diag.X(), center_diag.Y(), center_diag.Z());
	Edge edge_diag1 = m_mesh->newEdge(node_opposit, nodeQuad_split);
	Edge edge_diag0 = m_mesh->newEdge(nodeQuad_split, Anode);
	Edge edge_diag2 = m_mesh->newEdge(nodeQuad_split, node_new);

	nodeQuad_split.add<Edge>(edge_diag1);
	Anode.add<Edge>(edge_diag1);
	nodeQuad_split.add<Edge>(edge_diag0);
	Anode.add<Edge>(edge_diag0);
	nodeQuad_split.add<Edge>(edge_diag2);
	node_new.add<Edge>(edge_diag2);

	quadBig.remove<Node>(node_max);
	node_max.remove<Face>(quadBig);
	quadBig.remove<Edge>(edge_side);
	edge_side.remove<Face>(quadBig);
	quadBig.remove<Edge>(Aedge_max);
	Aedge_max.remove<Face>(quadBig);


	quadBig.add<Node>(nodeQuad_split);
	nodeQuad_split.add<Face>(quadBig);
	edge_diag0.add<Face>(quadBig);

	edge_diag1.add<Face>(quadBig);
	quadBig.add<Edge>(edge_diag1);


	std::vector<Face> listInside;
	std::vector<Edge> listIntersect = interestionS(node_new, node_min);
	Edge edgeRecovery = recoveryEdge(listIntersect, node_new, node_min);


	Face quad_new0 = createQuad(Anode, Aedge_min, edge_diag0, edgeRecovery, edge_diag2);
	std::cout<<quad_new0<<std::endl;
	Face quad_new1 = createQuad(node_new, Aedge_new1, edge_diag2, edge_side, edge_diag1);
	std::cout<<quad_new1<<std::endl;

	removeEdge(Aedge_max);

	return edgeRecovery;

}


/*----------------------------------------------------------------------------*/

void Quadfront::seam(Node& Anode, std::vector<Node>& listNode){
	double e1 = M_PI/6;
	double e2 = M_PI/4;

	Edge edge0 = std::get<1>(m_nodeFront[Anode]);
	Node node0 = edge0.getOppositeNode(Anode);
	Edge edge1 = std::get<2>(m_nodeFront[Anode]);
	Node node1 = edge1.getOppositeNode(Anode);
	double alpha = angle(Anode);

	if (edge0.get<Face>().size() == 1 or edge1.get<Face>().size() == 1) {
	}

	if ((alpha<e1 and Anode.get<Edge>().size() >= 5) or alpha<e2){
		std::cout<<"closing seam"<<std::endl;

		Face triangle0 = edge0.get<Face>()[0].nbEdges() == 4 ? edge0.get<Face>()[1] : edge0.get<Face>()[0];
		Edge edge_inter = get_opposit2Node(Anode, triangle0);
		Face triangle1 = edge_inter.get<Face>()[0] != triangle0 ? edge_inter.get<Face>()[0] : edge_inter.get<Face>()[1];

		Edge next_Aedge1 = std::get<1>(m_nodeFront[node1]) == edge1 ? std::get<2>(m_nodeFront[node1]) : std::get<1>(m_nodeFront[node1]);
		Face triangle_next_Aedge1 = next_Aedge1.get<Face>()[0].nbEdges() == 4 ? next_Aedge1.get<Face>()[1] : next_Aedge1.get<Face>()[0];
		Edge next_Aedge0 = std::get<1>(m_nodeFront[node0]) == edge0 ? std::get<2>(m_nodeFront[node0]) : std::get<1>(m_nodeFront[node0]);
		Face triangle_next_Aedge0 = next_Aedge0.get<Face>()[0].nbEdges() == 4 ? next_Aedge0.get<Face>()[1] : next_Aedge0.get<Face>()[0];

		Node node_top = get_opposit2Edge(edge_inter, triangle1);

		if(node_top == next_Aedge0.getOppositeNode(node0)){
			Edge split_inter = get_opposit2Node(node0, triangle1);
			Node new_node = m_mesh->newNode(split_inter.center());
			Face face_opposit = split_inter.get<Face>()[0] == triangle1 ? split_inter.get<Face>()[1] : split_inter.get<Face>()[0];
			split_operation(new_node, node0, triangle1, face_opposit);
		}
		else if (node_top == next_Aedge1.getOppositeNode(node1)){
			Edge split_inter = get_opposit2Node(node1, triangle1);
			Node new_node = m_mesh->newNode(split_inter.center());
			Face face_opposit = split_inter.get<Face>()[0] == triangle1 ? split_inter.get<Face>()[1] : split_inter.get<Face>()[0];
			split_operation(new_node, node1, triangle1, face_opposit);
		}

		triangle_next_Aedge1 = next_Aedge1.get<Face>()[0].nbEdges() == 4 ? next_Aedge1.get<Face>()[1] : next_Aedge1.get<Face>()[0];
		triangle_next_Aedge0 = next_Aedge0.get<Face>()[0].nbEdges() == 4 ? next_Aedge0.get<Face>()[1] : next_Aedge0.get<Face>()[0];
		triangle1 = edge_inter.get<Face>()[0] != triangle0 ? edge_inter.get<Face>()[0] : edge_inter.get<Face>()[1];
		node_top = get_opposit2Edge(edge_inter, triangle1);

		if(get_opposit2Edge(next_Aedge0, triangle_next_Aedge0) != node_top){
			Node next_node = next_Aedge0.getOppositeNode(node0);
			Edge swap_inter = get_opposit2Node(next_node, triangle_next_Aedge0);
			Face swap_face = swap_inter.get<Face>()[0] != triangle_next_Aedge0 ? swap_inter.get<Face>()[0] : swap_inter.get<Face>()[1];
			swap(next_node, triangle_next_Aedge0, swap_face);
		}

		if(get_opposit2Edge(next_Aedge1, triangle_next_Aedge1) != node_top){
			Node next_node = next_Aedge1.getOppositeNode(node1);
			Edge swap_inter = get_opposit2Node(next_node, triangle_next_Aedge1);
			Face swap_face = swap_inter.get<Face>()[0] != triangle_next_Aedge1 ? swap_inter.get<Face>()[0] : swap_inter.get<Face>()[1];
			swap(next_node, triangle_next_Aedge1, swap_face);
		}

		std::vector<Face> listInside;
		std::vector<Edge> listIntersect = interestionS(node0, node1);
		Edge edgeRecovery = recoveryEdge(listIntersect, node0, node1);
		math::Point midPoint = edgeRecovery.center();

		Edge edge_top;
		Edge edge_side;

		for (auto& e : node_top.get<Edge>()){
			if (e.getOppositeNode(node_top) == node0){
				edge_side = e;
			}
			else if (e.getOppositeNode(node_top) == node1){
				edge_top = e;
			}
		}

		Face quad_eph = createQuad(Anode, edge0, edge_side, edge1, edge_top);
		removeFace(quad_eph);

		node0.X() = midPoint.X();
		node0.Y() = midPoint.Y();


		Edge new_edge = closingSeam(edge0, edge1, edge_side, edge_top);

		listNode.push_back(Anode);
		listNode.push_back(node0);
		listNode.push_back(edge1.getOppositeNode(Anode));
	}

	else if(std::max(edge0.length(), edge1.length())/std::min(edge0.length(), edge1.length()) > 2.5){
		Edge edge_max = edge0.length() > edge1.length() ? edge0 : edge1;
		Edge edge_min = edge_max == edge0 ? edge1 : edge0 ;
		Edge edge_recovery = transitionSeam(Anode, edge_max, edge_min);

		Face quad_recovery = edge_recovery.get<Face>()[0].nbEdges() == 4 ? edge_recovery.get<Face>()[0] : edge_recovery.get<Face>()[1];
		Edge edge_baseRecovery = get_oppositEdge_Quad(edge_recovery, quad_recovery);
		update_Boundary2(edge_baseRecovery, edge_recovery);

		Node node0 = edge_recovery.get<Node>()[0];
		Node node1 = edge_recovery.get<Node>()[1];
		Edge edge0 = std::get<1>(m_nodeFront[node0]) == edge_recovery ? std::get<2>(m_nodeFront[node0]) : std::get<1>(m_nodeFront[node0]);
		Edge edge1 = std::get<1>(m_nodeFront[node1]) == edge_recovery ? std::get<2>(m_nodeFront[node1]) : std::get<1>(m_nodeFront[node1]);

		Face quad0 = edge0.get<Face>()[0].nbEdges() == 4 ? edge0.get<Face>()[0] : edge0.get<Face>()[1];
		Edge edge_base0 = get_oppositEdge_Quad(edge0, quad0);
		update_Boundary2(edge_base0, edge0);

		Face quad1 = edge1.get<Face>()[0].nbEdges() == 4 ? edge1.get<Face>()[0] : edge1.get<Face>()[1];
		Edge edge_base1 = get_oppositEdge_Quad(edge1, quad1);
		listNode.push_back(Anode);
		update_Boundary2(edge_base1, edge1);
	}
	else if (std::max(edge0.length(), edge1.length())/std::min(edge0.length(), edge1.length()) > 2.5 and alpha>e1 and alpha>e2){
		Edge edge_max = edge0.length() > edge1.length() ? edge0 : edge1;
		Face quad_max = edge_max.get<Face>()[0].nbEdges() == 4 ? edge_max.get<Face>()[0] : edge_max.get<Face>()[1];
		Edge edge_min = edge_max == edge0 ? edge1 : edge0 ;
		Face quad_min = edge_min.get<Face>()[0].nbEdges() == 4 ? edge_min.get<Face>()[0] : edge_min.get<Face>()[1];

		Edge edge_recovery = transitionSplit(Anode, edge_max, edge_min);


		Face quad_recovery = edge_recovery.get<Face>()[0].nbEdges() == 4 ? edge_recovery.get<Face>()[0] : edge_recovery.get<Face>()[1];
		Edge edge_baseRecovery = get_oppositEdge_Quad(edge_recovery, quad_recovery);
		update_Boundary2(edge_baseRecovery, edge_recovery);

		Node node0 = edge_recovery.get<Node>()[0] == edge_min.getOppositeNode(Anode) ? edge_recovery.get<Node>()[1] : edge_recovery.get<Node>()[0];
		Edge edge0 = std::get<1>(m_nodeFront[node0]) == edge_recovery ? std::get<2>(m_nodeFront[node0]) : std::get<1>(m_nodeFront[node0]);
		Face face0 = edge0.get<Face>()[0].nbEdges() == 4 ? edge0.get<Face>()[0] : edge0.get<Face>()[1];
		Edge edge_base0 = get_oppositEdge_Quad(edge0, face0);
		listNode.push_back(Anode);
		update_Boundary2(edge_base0, edge0);
	}

}

/*----------------------------------------------------------------------------*/

Face Quadfront::qmorphAlgo(Edge& Aedge){

	std::cout<<"\n======= Edge "<<Aedge.id()<<"========\n"<<std::endl;
	std::cout<<Aedge.get<Node>()[0]<<std::endl;
	std::cout<<Aedge.get<Node>()[1]<<std::endl;

	Node node0 =Aedge.get<Node>()[0];
	Node node1 =Aedge.get<Node>()[1];
	Edge sideEdge0;
	Node nodeRecovery0;

	std::cout<<"\n======= SideEdge "<<node0<<" ========\n"<<std::endl;

	std::cout<<"Status "<<node0<<" : "<<std::get<0>(m_nodeFront[node0])<<std::endl;
	bool test0 = false;
	if (std::get<0>(m_nodeFront[node0]) == 1){
		sideEdge0 = (std::get<1>(m_nodeFront[node0]) == Aedge) ? std::get<2>(m_nodeFront[node0]) : std::get<1>(m_nodeFront[node0]);
	}
	else{
		sideEdge0 = sideEdge(node0, Aedge);
	}

	std::cout<<"\nsideEdge0 "<<sideEdge0.id()<<" : "<<sideEdge0.get<Node>()[0]<<" ; "<<sideEdge0.get<Node>()[1]<<std::endl;
	std::cout<<"nodeRecovery0 : "<<sideEdge0.getOppositeNode(node0)<<std::endl;

	std::cout<<"\n======= SideEdge "<<node1<<"  ========\n"<<std::endl;
	Edge sideEdge1;
	Node nodeRecovery1;

	std::cout<<"Status "<<node1<<" : "<<std::get<0>(m_nodeFront[node1])<<std::endl;

	bool test1 = false;
	if (std::get<0>(m_nodeFront[node1]) == 1){
		sideEdge1 = (std::get<1>(m_nodeFront[node1]) == Aedge) ? std::get<2>(m_nodeFront[node1]) : std::get<1>(m_nodeFront[node1]);
	}
	else{
		sideEdge1 = sideEdge(node1, Aedge);
	}

	std::cout<<"\nsideEdge1 "<<sideEdge1.id()<<" : "<<sideEdge1.get<Node>()[0]<<" ; "<<sideEdge1.get<Node>()[1]<<std::endl;
	std::cout<<"nodeRecovery1 : "<<sideEdge1.getOppositeNode(node1)<<std::endl;



	// CAS PARTICULIER nodeRecovery0 = nodeRecovery1
	if (sideEdge1.getOppositeNode(node1) == sideEdge0.getOppositeNode(node0)){
		Face face_com;
		for (auto f : Aedge.get<Face>()){
			if (f.nbEdges() == 3){
				face_com = f;
			}
		}

		if (std::get<0>(m_nodeFront[node0]) == 0){
			std::cout<<"\n======= Modification SideEdge0 ========\n"<<std::endl;
			Face face_split = (sideEdge0.get<Face>()[0] == face_com) ? sideEdge0.get<Face>()[1] : sideEdge0.get<Face>()[0];
			Edge edge_split_inter = get_opposit2Node(node0, face_split);
			Face face_split_opposit = (edge_split_inter.get<Face>()[0] == face_split) ? edge_split_inter.get<Face>()[1] : edge_split_inter.get<Face>()[0];
			Node node_split_opposit = get_opposit2Edge(edge_split_inter, face_split_opposit);

			std::map<Node, std::pair<Edge, Edge>> mapFace_split = get_edgeInFace(face_split);
			std::map<Node, std::pair<Edge, Edge>> mapFace_split_opposit = get_edgeInFace(face_split_opposit);

			Edge edge0 = std::get<0>(mapFace_split[node0]);
			Edge edge1 = std::get<1>(mapFace_split[node0]);

			math::Vector3d bissectrice_split = get_bissectrice(node0, edge0, edge1);

			math::Vector3d vNode = intersectionVec2Edge(bissectrice_split, node0, edge_split_inter);
			Node newNode = m_mesh->newNode(vNode.X(), vNode.Y(), vNode.Z());

			sideEdge0 = split_operation(newNode, node0, face_split, face_split_opposit);

			std::cout<<"\nsideEdge0 "<<sideEdge0.id()<<" : "<<sideEdge0.get<Node>()[0]<<" ; "<<sideEdge0.get<Node>()[1]<<std::endl;
			std::cout<<"nodeRecovery0 : "<<sideEdge0.getOppositeNode(node0)<<std::endl;

		}
		else if (std::get<0>(m_nodeFront[node1]) == 0){
			std::cout<<"\n======= Modification SideEdge1 ========\n"<<std::endl;
			Face face_split = (sideEdge1.get<Face>()[0] == face_com) ? sideEdge1.get<Face>()[1] : sideEdge1.get<Face>()[0];
			Edge edge_split_inter = get_opposit2Node(node1, face_split);
			Face face_split_opposit = (edge_split_inter.get<Face>()[0] == face_split) ? edge_split_inter.get<Face>()[1] : edge_split_inter.get<Face>()[0];
			Node node_split_opposit = get_opposit2Edge(edge_split_inter, face_split_opposit);

			std::map<Node, std::pair<Edge, Edge>> mapFace_split = get_edgeInFace(face_split);
			std::map<Node, std::pair<Edge, Edge>> mapFace_split_opposit = get_edgeInFace(face_split_opposit);

			Edge edge0 = std::get<0>(mapFace_split[node1]);
			Edge edge1 = std::get<1>(mapFace_split[node1]);

			math::Vector3d bissectrice_split = get_bissectrice(node1, edge0 , edge1);

			math::Vector3d vNode = intersectionVec2Edge(bissectrice_split, node1, edge_split_inter);
			Node newNode = m_mesh->newNode(vNode.X(), vNode.Y(), vNode.Z());

			sideEdge1 = split_operation(newNode, node1, face_split, face_split_opposit);

			std::cout<<"\nsideEdge1 "<<sideEdge1.id()<<" : "<<sideEdge1.get<Node>()[0]<<" ; "<<sideEdge1.get<Node>()[1]<<std::endl;
			std::cout<<"nodeRecovery1 : "<<sideEdge1.getOppositeNode(node1)<<std::endl;

		}
	}

	nodeRecovery0 = sideEdge0.getOppositeNode(node0);
	nodeRecovery1 = sideEdge1.getOppositeNode(node1);


	math::Vector3d vecAB({nodeRecovery1.X()-nodeRecovery0.X(),
	                      nodeRecovery1.Y()-nodeRecovery0.Y(),
	                      nodeRecovery1.Z()-nodeRecovery0.Z()});


	std::cout<<"\n ==== UpSide ====\n "<<std::endl;

	Edge edgeRecovery;

	std::vector<Edge> listIntersect = interestionS(nodeRecovery0, nodeRecovery1);

	std::cout<<"listIntersect "<<listIntersect.size()<<std::endl;

	edgeRecovery = recoveryEdge(listIntersect, nodeRecovery0, nodeRecovery1);


	std::cout<<"edgeRecovery : "<<edgeRecovery.get<Node>()[0]<<" ; "<<edgeRecovery.get<Node>()[1]<<std::endl;
	Face quad = createQuad(node0, Aedge, sideEdge0, sideEdge1, edgeRecovery);

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.setDataOptions(gmds::N | gmds::F | gmds::E);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link2.vtk");

	update_Boundary2(Aedge, edgeRecovery);

	return quad;
}

/*----------------------------------------------------------------------------*/

std::vector<double> Quadfront::statFront(std::vector<Edge> listEdge){
	std::vector<double> stat;
	double max;
	double min;
	double av = 0;
	if (listEdge[0].length() >= listEdge[1].length()){
		max = listEdge[0].length();
		min = listEdge[1].length();
	}
	else{
		max = listEdge[1].length();
		min = listEdge[0].length();
	}
	av += (max+min) + listEdge.size();

	for (auto& i : listEdge){
		if (i != listEdge[0] and i != listEdge[1]){
			if (i.length() > max){
				max = i.length();
			}
			else if (i.length() < min){
				min = i.length();
			}
			av += i.length()/listEdge.size();
		}
	}
	stat.push_back(av);
	std::cout<<" average length : "<<av<<std::endl;
	stat.push_back(max/min);
	std::cout<<" max/min : "<<max/min<<std::endl;
	return stat;
}

/*----------------------------------------------------------------------------*/

void Quadfront::qmorphSmoothing(Node& Anode){


	double length;
	math::Vector3d vec;

	Edge edge_up0 = std::get<1>(m_nodeFront[Anode]);
	Edge edge_up1 = std::get<2>(m_nodeFront[Anode]);

	Face quad0 = edge_up0.get<Face>()[0].nbNodes() == 4 ? edge_up0.get<Face>()[0] : edge_up0.get<Face>()[1];
	Face quad1 = edge_up1.get<Face>()[0].nbNodes() == 4 ? edge_up1.get<Face>()[0] : edge_up1.get<Face>()[1];
	std::map<Node, std::pair<Edge, Edge>> mapFace0 = get_edgeInFace(quad0);
	std::map<Node, std::pair<Edge, Edge>> mapFace1 = get_edgeInFace(quad1);

	Node opposit_node0 = get_oppositNode_Quad(Anode, quad0);
	Node opposit_node1 = get_oppositNode_Quad(Anode, quad1);


	length = std::get<0>(mapFace0[opposit_node0]).length() +
	         std::get<1>(mapFace0[opposit_node0]).length() +
	         std::get<0>(mapFace1[opposit_node1]).length() +
	         std::get<1>(mapFace1[opposit_node1]).length();

	int numb=0;

	for (auto& e : Anode.get<Edge>()){
		if (e != std::get<0>(mapFace0[Anode]) and e != std::get<1>(mapFace0[Anode]) and e != std::get<0>(mapFace1[Anode]) and e != std::get<1>(mapFace1[Anode])){
			length += e.length();
			numb+=1;
		}
	}

	length = length/(4+numb);

	Edge sideEdge = std::get<0>(mapFace0[Anode]) != edge_up0 ? std::get<0>(mapFace0[Anode]) : std::get<1>(mapFace0[Anode]);
	Node node_interior = sideEdge.getOppositeNode(Anode);

	math::Vector3d vec_SideEdge = node2Vector(sideEdge, node_interior).normalize();
	math::Vector3d vecEdge = length*vec_SideEdge;

	std::cout<<"Ancien Anode "<<Anode<<std::endl;


	Anode.X() = node_interior.X() + vecEdge.X();
	Anode.Y() = node_interior.Y() + vecEdge.Y();
	Anode.Z() = node_interior.Z() + vecEdge.Z();

	std::cout<<"Nouveau Anode "<<Anode<<std::endl;
}

/*----------------------------------------------------------------------------*/

void Quadfront::blSmoothing(Node& Anode, double average, double tr){
	double length;
	math::Point vec = {0,  0, 0};
	int i = 0;
	Edge edge_up;
	Edge edge_side;
	Node node_interior;
	for (auto& f : Anode.get<Face>()){
		if (f.nbEdges() == 4){
			std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(f);
			if (std::get<0>(mapFace[Anode]) == std::get<2>(m_nodeFront[Anode]) or std::get<0>(mapFace[Anode]) == std::get<1>(m_nodeFront[Anode])){
				edge_up = std::get<0>(mapFace[Anode]);
				edge_side = std::get<1>(mapFace[Anode]);
			}
			else{
				edge_up = std::get<1>(mapFace[Anode]);
				edge_side = std::get<0>(mapFace[Anode]);
			}
			node_interior = edge_side.getOppositeNode(Anode);
			length = edge_side.length();
			vec.X() += edge_up.getOppositeNode(Anode).X() +
						  node_interior.X() - get_oppositNode_Quad(Anode, f).X();
			vec.Y() += edge_up.getOppositeNode(Anode).Y() +
						  node_interior.Y() - get_oppositNode_Quad(Anode, f).Y();
			i++;
		}
	}
	vec = {vec.X()/i, vec.Y()/i, vec.Z()/i};
	math::Point deltaA = {vec.X() - Anode.X(), vec.Y() - Anode.Y(), vec.Z() - Anode.Z()};
	Anode.X() = vec.X();
	Anode.Y() = vec.Y();

	math::Point deltaB;
	deltaB.X() = node2Vector(edge_side, Anode).X() + (deltaA.X() + node2Vector(edge_side, node_interior).X())*length/edge_side.length();
	deltaB.Y() = node2Vector(edge_side, Anode).Y() + (deltaA.Y() + node2Vector(edge_side, node_interior).Y())*length/edge_side.length();
	Anode.X() += deltaB.X();
	Anode.Y() += deltaB.Y();

	Face quad0 = edge_side.get<Face>()[0];
	Node node0 = get_oppositNode_Quad(node_interior, quad0);
	Face quad1 = edge_side.get<Face>()[1];
	Node node1 = get_oppositNode_Quad(node_interior, quad1);

	math::Vector3d pb_quad0 = {node0.X() - node_interior.X(), node0.Y() - node_interior.Y(), node0.Z() - node_interior.Z()};
	math::Vector3d pb_quad1 = {node1.X() - node_interior.X(), node1.Y() - node_interior.Y(), node1.Z() - node_interior.Z()};
	math::Vector3d pb1 = pb_quad0 + pb_quad1;

	math::Vector3d pb2 = pb1.normalize()+ node2Vector(edge_side, node_interior).normalize();

	Edge edge_q = m_mesh->newEdge(node0, node1);
	math::Vector3d vecQ = intersectionVec2Edge(pb2.normalize(), node_interior, edge_q);
	m_mesh->deleteEdge(edge_q);
	vecQ = {vecQ.X()-node_interior.X(), vecQ.Y()-node_interior.Y(), vecQ.Z()-node_interior.Z()};

	if(average > 20 and tr > 2.5){
		qmorphSmoothing(Anode);
	}

	double length_pb2 = length > vecQ.norm() ? (vecQ.norm() + length)/2 : length;

	pb2 = pb2.normalize()*length_pb2;

	Anode.X() = node_interior.X() + pb2.X();
	Anode.Y() = node_interior.Y() + pb2.Y();

}

/*----------------------------------------------------------------------------*/

void Quadfront::interiorSmoothing(Node& Anode){
	math::Vector3d vec;
	for (auto &e : Anode.get<Edge>()) {
		vec = vec + node2Vector(e, Anode)/ Anode.get<Edge>().size();
	}
	Anode.X() += vec.X() ;
	Anode.Y() += vec.Y() ;
	Anode.Z() += vec.Z() ;
}

/*----------------------------------------------------------------------------*/


void Quadfront::triangleSmoothing(Edge Aedge){

	for (auto &node : Aedge.get<Node>()) {
		for (auto &e0 : node.get<Edge>()) {
			std::vector<math::Point> vec_point;
			Node opposit_node = e0.getOppositeNode(node);

			if (e0.get<Face>().size() == 2) {
				if (e0.get<Face>()[0].nbEdges() != 4 and e0.get<Face>()[1].nbEdges() != 4 and m_nodeFront.find(opposit_node) == m_nodeFront.end()) {
					for (auto &ee0 : opposit_node.get<Edge>()) {
						math::Point node_for = ee0.getOppositeNode(opposit_node).point();
						if (std::find(vec_point.begin(), vec_point.end(), node_for) == vec_point.end()) {
							vec_point.push_back(node_for);
						}
					}
					math::Point center = math::Point::massCenter(vec_point);
					opposit_node.X() = center.X();
					opposit_node.Y() = center.Y();
					opposit_node.Z() = center.Z();
				}
			}
		}
	}

}

/*----------------------------------------------------------------------------*/

void Quadfront::initCorners(Node& Anode, std::vector<Face>& listQuad){
	math::Vector3d bis_vec = -get_bissectrice(Anode, std::get<0>(m_nodeBoundary[Anode]), std::get<1>(m_nodeBoundary[Anode])).normalize();
	Face bestFace = Anode.get<Face>()[0];
	double dot = 0;
	for(auto& f : Anode.get<Face>()){
		std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(f);
		math::Vector3d vec = get_bissectrice(Anode, std::get<0>(mapFace[Anode]), std::get<1>(mapFace[Anode])).normalize();
		if (bis_vec.dot(vec) > dot){
			bestFace = f;
			dot = bis_vec.dot(vec);
		}
	}

	std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(bestFace);
	Edge edge_base = std::get<0>(mapFace[Anode]);
	Node node_base = edge_base.getOppositeNode(Anode);
	Edge edge_inter = get_opposit2Node(Anode, bestFace);
	Face face_opposit = edge_inter.get<Face>()[0] == bestFace ? edge_inter.get<Face>()[1] : edge_inter.get<Face>()[0];
	Node node_opposit = get_opposit2Edge(edge_inter, face_opposit);
	std::map<Node, std::pair<Edge, Edge>> mapFace_opposit = get_edgeInFace(face_opposit);
	Edge edge_side;
	Edge edge_top;

	std::cout<<std::get<0>(mapFace_opposit[node_opposit]).getOppositeNode(node_opposit)<<std::endl;
	std::cout<<edge_base.getOppositeNode(Anode)<<std::endl;

	if (std::get<0>(mapFace_opposit[node_opposit]).getOppositeNode(node_opposit) == node_base){
		edge_side = std::get<0>(mapFace_opposit[node_opposit]);
		edge_top = std::get<1>(mapFace_opposit[node_opposit]);
	}
	else{
		edge_side = std::get<1>(mapFace_opposit[node_opposit]);
		edge_top = std::get<0>(mapFace_opposit[node_opposit]);
	}

	Face quad = createQuad(Anode, edge_base, edge_side, std::get<1>(mapFace[Anode]), edge_top);
	removeEdge(edge_inter);
	update_Boundary2(edge_base, edge_top);
	listQuad.push_back(quad);

	if(edge_base.get<Face>()[0].nbEdges() == 3 or edge_base.get<Face>()[1].nbEdges() == 3) {
		Face face0 = edge_base.get<Face>()[0] == quad ? edge_base.get<Face>()[1] : edge_base.get<Face>()[0];
		Edge edge_top0 = get_opposit2Node(Anode, face0);
		Edge edge_inter0 = get_opposit2Node(node_base, face0);
		Face face0_opposit = edge_inter0.get<Face>()[0] == face0 ? edge_inter0.get<Face>()[1] : edge_inter0.get<Face>()[0];
		Edge edge_side0 = get_opposit2Node(Anode, face0_opposit);
		Edge edge_base0 =
		   std::get<0>(m_nodeBoundary[Anode]).get<Face>()[0] == face0_opposit ? std::get<0>(m_nodeBoundary[Anode]) : std::get<1>(m_nodeBoundary[Anode]);

		Face quad0 = createQuad(Anode, edge_base0, edge_side0, edge_base, edge_top0);
		update_Boundary2(edge_base0, edge_top0);
		listQuad.push_back(quad0);
	}

	if(std::get<1>(mapFace[Anode]).get<Face>()[0].nbEdges() == 3 or std::get<1>(mapFace[Anode]).get<Face>()[1].nbEdges() == 3) {
		Face face1 = std::get<1>(mapFace[Anode]).get<Face>()[1] == quad ? std::get<1>(mapFace[Anode]).get<Face>()[0] : std::get<1>(mapFace[Anode]).get<Face>()[1];
		Edge edge_top1 = get_opposit2Node(Anode, face1);
		Node node_side = std::get<1>(mapFace[Anode]).getOppositeNode(Anode);
		Edge edge_inter1 = get_opposit2Node(node_side, face1);
		Face face1_opposit = edge_inter1.get<Face>()[0] == face1 ? edge_inter1.get<Face>()[1] : edge_inter1.get<Face>()[0];
		Edge edge_side1 = get_opposit2Node(Anode, face1_opposit);
		Edge edge_base1 =
		   std::get<0>(m_nodeBoundary[Anode]).get<Face>()[0] == face1_opposit ? std::get<0>(m_nodeBoundary[Anode]) : std::get<1>(m_nodeBoundary[Anode]);

		Face quad1 = createQuad(Anode, edge_base1, edge_side1, std::get<1>(mapFace[Anode]), edge_top1);
		update_Boundary2(edge_base1, edge_top1);
		listQuad.push_back(quad1);
	}
}


/*----------------------------------------------------------------------------*/

void Quadfront::petitAngle(Face &Aface, std::vector<Face> &listFace, std::vector<Edge> &listEdge){
	if (std::find(listFace.begin(), listFace.end(), Aface) != listFace.end() or Aface.nbEdges() == 4) {
		return;
	}
	listFace.push_back(Aface);
	bool test = true;
	Node node;
	for(auto& n : Aface.get<Node>()){
		if (m_nodeBoundary.find(n) == m_nodeBoundary.end()){
			test = false;
			node = n;
			break;
		}
	}
	if (test) {
		for (auto& e : Aface.get<Edge>()){
			if (std::get<1>(m_nodeBoundary[e.get<Node>()[0]]) != e and std::get<0>(m_nodeBoundary[e.get<Node>()[0]]) != e and
			    std::get<1>(m_nodeBoundary[e.get<Node>()[1]]) != e and std::get<0>(m_nodeBoundary[e.get<Node>()[1]]) != e){
				if (std::find(listEdge.begin(), listEdge.end(), e) == listEdge.end()){
					listEdge.push_back(e);
					Face face = e.get<Face>()[0] == Aface ? e.get<Face>()[1] : e.get<Face>()[0];
					petitAngle(face, listFace, listEdge);
					return;
				}
			}
		}
	}
	else{
		for(auto& ea : node.get<Edge>()){
			if (m_nodeFront.find(ea.getOppositeNode(node)) == m_nodeFront.end()){
				return;
			}
		}
		for(auto& fa : node.get<Face>()){
			if (fa != Aface){
				Edge edge_opposit = get_opposit2Node(node, fa);
				if (std::get<1>(m_nodeBoundary[edge_opposit.get<Node>()[0]]) != edge_opposit and
				    std::get<0>(m_nodeBoundary[edge_opposit.get<Node>()[0]]) != edge_opposit and
				    std::get<1>(m_nodeBoundary[edge_opposit.get<Node>()[1]]) != edge_opposit and
				    std::get<0>(m_nodeBoundary[edge_opposit.get<Node>()[1]]) != edge_opposit){
					listEdge.push_back(edge_opposit);
					Face face = edge_opposit.get<Face>()[0] == Aface ? edge_opposit.get<Face>()[1] : edge_opposit.get<Face>()[0];
					if (std::find(listFace.begin(), listFace.end(), Aface) == listFace.end()) {
						petitAngle(face, listFace, listEdge);
					}
					return;
				}
			}
		}
		return;
	}
}

/*----------------------------------------------------------------------------*/

void Quadfront::initRevers(gmds::Node &Anode, std::vector<Face>& listQuad){
	math::Vector3d bissec = (-node2Vector(std::get<0>(m_nodeBoundary[Anode]), Anode).normalize() - node2Vector(std::get<1>(m_nodeBoundary[Anode]), Anode).normalize()).normalize();
	Edge edge_max = Anode.get<Edge>()[0];
	double dot_max = node2Vector(edge_max, Anode).normalize().dot(bissec);
	for(auto& e : Anode.get<Edge>()){
		math::Vector3d vec = node2Vector(e, Anode).normalize();
		if (vec.dot(bissec) > dot_max){
			dot_max = vec.dot(bissec);
			edge_max = e;
		}
	}

	Face face0 = edge_max.get<Face>()[0];
	Face face1 = edge_max.get<Face>()[1];
	Edge edge_inter0 = get_opposit2Node(Anode, face0);
	Face face0_opposit = edge_inter0.get<Face>()[0] == face0 ? edge_inter0.get<Face>()[1] : edge_inter0.get<Face>()[0];
	Node node0_opposit = get_opposit2Edge(edge_inter0, face0_opposit);
	std::map<Node, std::pair<Edge, Edge>> mapFace0 = get_edgeInFace(face0);
	std::map<Node, std::pair<Edge, Edge>> mapFace0_opposit = get_edgeInFace(face0_opposit);

	Edge edge0 = std::get<0>(mapFace0[Anode]) == edge_max ? std::get<1>(mapFace0[Anode]) : std::get<0>(mapFace0[Anode]);
	Node node_max = edge_max.getOppositeNode(Anode);
	Face quad0;

	if(std::get<0>(mapFace0_opposit[node0_opposit]).getOppositeNode(node0_opposit) == node_max) {
		quad0 = createQuad(Anode, edge_max, edge0, std::get<0>(mapFace0_opposit[node0_opposit]), std::get<1>(mapFace0_opposit[node0_opposit]));
	}
	else{
		quad0 = createQuad(Anode, edge_max, edge0, std::get<1>(mapFace0_opposit[node0_opposit]), std::get<0>(mapFace0_opposit[node0_opposit]));
	}

	Edge edge_top0 = get_oppositEdge_Quad(edge_max, quad0);
	update_Boundary2(edge_max, edge_top0);
	listQuad.push_back(quad0);

	Edge edge_inter1 = get_opposit2Node(Anode, face1);
	Face face1_opposit = edge_inter1.get<Face>()[0] == face1 ? edge_inter1.get<Face>()[1] : edge_inter1.get<Face>()[0];
	Node node1_opposit = get_opposit2Edge(edge_inter1, face1_opposit);
	std::map<Node, std::pair<Edge, Edge>> mapFace1 = get_edgeInFace(face1);
	std::map<Node, std::pair<Edge, Edge>> mapFace1_opposit = get_edgeInFace(face1_opposit);

	Edge edge1 = std::get<0>(mapFace1[Anode]) == edge_max ? std::get<1>(mapFace1[Anode]) : std::get<0>(mapFace1[Anode]);

	Face quad1;

	if(std::get<0>(mapFace1_opposit[node1_opposit]).getOppositeNode(node1_opposit) == node_max) {
		quad1 = createQuad(Anode, edge_max, edge1, std::get<0>(mapFace1_opposit[node1_opposit]), std::get<1>(mapFace1_opposit[node1_opposit]));
	}
	else{
		quad1 = createQuad(Anode, edge_max, edge1, std::get<1>(mapFace1_opposit[node1_opposit]), std::get<0>(mapFace1_opposit[node1_opposit]));
	}

	Edge edge_top1 = get_oppositEdge_Quad(edge_max, quad1);
	update_Boundary2(edge_max, edge_top1);
	listQuad.push_back(quad1);

	if(edge1.get<Face>()[0].nbEdges() == 3 or edge1.get<Face>()[1].nbEdges() == 3) {
		Face face2 = edge1.get<Face>()[1] == quad1 ? edge1.get<Face>()[0] : edge1.get<Face>()[1];
		Edge edge_top2 = get_opposit2Node(Anode, face2);
		Node node_side2 = edge1.getOppositeNode(Anode);
		Edge edge_inter2 = get_opposit2Node(node_side2, face2);
		Face face2_opposit = edge_inter2.get<Face>()[0] == face2 ? edge_inter2.get<Face>()[1] : edge_inter2.get<Face>()[0];
		Edge edge_side2 = get_opposit2Node(Anode, face2_opposit);
		Edge edge_base2 =
		   std::get<0>(m_nodeBoundary[Anode]).get<Face>()[0] == face2_opposit ? std::get<0>(m_nodeBoundary[Anode]) : std::get<1>(m_nodeBoundary[Anode]);

		Face quad2 = createQuad(Anode, edge_base2, edge1, edge_side2, edge_top2);
		update_Boundary2(edge_base2, edge_top2);
		listQuad.push_back(quad2);
	}

	if(edge0.get<Face>()[0].nbEdges() == 3 or edge0.get<Face>()[1].nbEdges() == 3) {
		Face face3 = edge0.get<Face>()[1] == quad0 ? edge0.get<Face>()[0] : edge0.get<Face>()[1];
		Edge edge_top3 = get_opposit2Node(Anode, face3);
		Node node_side3 = edge0.getOppositeNode(Anode);
		Edge edge_inter3 = get_opposit2Node(node_side3, face3);
		Face face3_opposit = edge_inter3.get<Face>()[0] == face3 ? edge_inter3.get<Face>()[1] : edge_inter3.get<Face>()[0];
		Edge edge_side3 = get_opposit2Node(Anode, face3_opposit);
		Edge edge_base3 =
		   std::get<0>(m_nodeBoundary[Anode]).get<Face>()[0] == face3_opposit ? std::get<0>(m_nodeBoundary[Anode]) : std::get<1>(m_nodeBoundary[Anode]);

		Face quad3 = createQuad(Anode, edge_base3, edge0, edge_side3, edge_top3);
		update_Boundary2(edge_base3, edge_top3);
		listQuad.push_back(quad3);
	}

}

/*----------------------------------------------------------------------------*/

void Quadfront::initEnds(Node& Anode, std::vector<Face>& listQuad){
	std::vector<Face> listFace;
	std::vector<Node> listNode;
	std::vector<Edge> listEdge;
	petitAngle(Anode.get<Face>()[0], listFace, listEdge);
	std::cout<<"listEdge : "<<listEdge.size()<<std::endl;
	if (listEdge.size()>=1) {
		for (auto i = 0; i != floor(listEdge.size()/2); i++) {
			Edge edge = listEdge[2*i+1];
			Node node_new = m_mesh->newNode(edge.center());
			Node Anode = get_opposit2Edge(edge, edge.get<Face>()[0]);
			Edge side_edge = split_operation(node_new, Anode, edge.get<Face>()[0], edge.get<Face>()[1]);
		}
		for (auto i = 1; i != floor(listEdge.size()/2); i++) {
			Edge edge = listEdge[2*i];
			Node node_new = m_mesh->newNode(edge.center());
			Node Anode = get_opposit2Edge(edge, edge.get<Face>()[0]);
			Edge side_edge = split_operation(node_new, Anode, edge.get<Face>()[0], edge.get<Face>()[1]);
		}
	}

	if (listEdge.size() % 2 == 1){
		Edge edge = listEdge.back();
		std::cout<<edge.get<Node>()[0]<<" ; "<<edge.get<Node>()[1]<<std::endl;
		Node node_new = m_mesh->newNode(edge.center());
		Node Anode = get_opposit2Edge(edge, edge.get<Face>()[0]);
		Edge side_edge = split_operation(node_new, Anode, edge.get<Face>()[0], edge.get<Face>()[1]);
	}

	Node node0 = std::get<0>(m_nodeBoundary[Anode]).getOppositeNode(Anode);
	Node node1 = std::get<1>(m_nodeBoundary[Anode]).getOppositeNode(Anode);
	Edge edge0 = std::get<0>(m_nodeBoundary[node0]) != std::get<0>(m_nodeBoundary[Anode]) and
	                   std::get<0>(m_nodeBoundary[node0]) != std::get<1>(m_nodeBoundary[Anode]) ?
	                std::get<0>(m_nodeBoundary[node0]) : std::get<1>(m_nodeBoundary[node0]);
	Edge edge1 = std::get<0>(m_nodeBoundary[node1]) != std::get<0>(m_nodeBoundary[Anode]) and
	                   std::get<0>(m_nodeBoundary[node1]) != std::get<1>(m_nodeBoundary[Anode]) ?
	                std::get<0>(m_nodeBoundary[node1]) : std::get<1>(m_nodeBoundary[node1]);

	Node node0_next = edge0.getOppositeNode(node0);
	Node node1_next = edge1.getOppositeNode(node1);

	Edge edge_top0 = get_opposit2Node(node0_next, edge0.get<Face>()[0]);
	Edge edge_top1 = get_opposit2Node(node1_next, edge1.get<Face>()[0]);

	Face quad = createQuad(Anode, std::get<0>(m_nodeBoundary[Anode]), std::get<1>(m_nodeBoundary[Anode]), edge_top0, edge_top1);

	update_Boundary2(std::get<0>(m_nodeBoundary[Anode]), edge_top1);
	listQuad.push_back(quad);
}

/*----------------------------------------------------------------------------*/

void Quadfront::initialiseQuad(gmds::Node &Anode, std::vector<Face>& listQuad){
	double sommeAngle = 0;
	for (auto &face : Anode.get<Face>()) {
		std::map<Node, std::pair<Edge, Edge>> node2edge = get_edgeInFace(face);
		sommeAngle += node2Vector(std::get<0>(node2edge[Anode]), Anode).normalize().angle(node2Vector(std::get<1>(node2edge[Anode]), Anode).normalize());

	}
	std::cout<<sommeAngle*360/(2*M_PI)<<std::endl;
	if(sommeAngle < M_PI/3) {
		initEnds(Anode, listQuad);
		/*
		Face quad = listQuad.back();
		std::map<Node, std::pair<Edge, Edge>> mapFace = get_edgeInFace(quad);
		Node node0 = std::get<1>(m_nodeBoundary[Anode]).getOppositeNode(Anode);
		if(std::get<1>(m_nodeFront[node0]) == std::get<1>(m_nodeFront[Anode])){
			std::get<1>(m_nodeFront[node0]) = std::get<1>(mapFace[node0]) == std::get<1>(m_nodeFront[node0]) ? std::get<0>(mapFace[node0]) : std::get<1>(mapFace[node0]);
		}
		else{
			std::get<2>(m_nodeFront[node0]) = std::get<1>(mapFace[node0]) == std::get<2>(m_nodeFront[node0]) ? std::get<0>(mapFace[node0]) : std::get<1>(mapFace[node0]);
		}
		std::get<0>(m_nodeFront[node0]) == 1;

		Node node1 = std::get<2>(m_nodeFront[Anode]).getOppositeNode(Anode);
		if(std::get<1>(m_nodeFront[node1]) == std::get<1>(m_nodeFront[Anode])){
			std::get<1>(m_nodeFront[node1]) = std::get<1>(mapFace[node1]) == std::get<1>(m_nodeFront[node1]) ? std::get<0>(mapFace[node1]) : std::get<1>(mapFace[node1]);
		}
		else{
			std::get<2>(m_nodeFront[node1]) = std::get<1>(mapFace[node1]) == std::get<2>(m_nodeFront[node1]) ? std::get<0>(mapFace[node1]) : std::get<1>(mapFace[node1]);
		}
		std::get<0>(m_nodeFront[node1]) == 1;
		 */
	}
	else if(sommeAngle > 5*M_PI/3) {
		initRevers(Anode, listQuad);
	}
	else if (sommeAngle < 5*M_PI/3 and sommeAngle > 4*M_PI/3){
		initCorners(Anode, listQuad);
	}
}

/*----------------------------------------------------------------------------*/
}     // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/