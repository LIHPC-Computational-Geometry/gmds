/*----------------------------------------------------------------------------*/
#include "gmds/quadfront/Quadfront.h"
#include "gmds/io/VTKReader.h"
#include <gmds/io/VTKWriter.h>
#include "gmds/ig/MeshDoctor.h"
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Triangle.h>
#include <cmath>

/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace quadfront {
/*----------------------------------------------------------------------------*/
Quadfront::Quadfront(const std::string& AFilename)
{

	m_mesh = new gmds::Mesh(gmds::MeshModel(gmds::DIM2 |  gmds::F | gmds::E | gmds::N |
	                                        gmds::N2E |gmds::N2F | gmds::E2N |
	                                        gmds::E2F | gmds::F2N | gmds::F2E));

	Variable<int>* nb_boundary = m_mesh->newVariable<int, GMDS_NODE>("nb_boundary");

	std::cout<<"============== Start Initialisation =============="<<std::endl;
	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	//vtkReader.setDataOptions(gmds::N | gmds::F);
	vtkReader.read(AFilename);

	gmds::MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();


	std::cout<<"Nombre de noeuds: " <<m_mesh->getNbNodes()<<std::endl;
	std::cout<<"Nombre d'arêtes: " <<m_mesh->getNbEdges()<<std::endl;
	std::cout<<"Nombre de faces: " <<m_mesh->getNbFaces()<<std::endl;

	std::vector<Edge> FrontEdges;


	for (auto& e : m_mesh->edges()){
		Edge edge = m_mesh->get<Edge>(e);
		if (edge.nbFaces() == 1){
			FrontEdges.push_back(edge);
		}
	}


	int i =1;
	while (!FrontEdges.empty()) {
		Edge e0 = FrontEdges.front();
		std::vector<Node> e0_nodes = e0.get<Node>();

		Node nEnd = e0_nodes[0];
		Node nStart = e0_nodes[1];
		m_nodeBoundary[nEnd] = std::make_tuple(0, e0, e0);
		m_edgeBoundary[e0] = i;

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
			std::set_intersection(vEdge.begin(),
			                      vEdge.end(),
			                      FrontEdges.begin(),
			                      FrontEdges.end(),
			                      std::back_inserter(intersection));
			if (size(intersection) == 0) { //Intersection is empty
				break;
			}
			else {
				m_nodeBoundary[nStart] = std::make_tuple(0, e0, intersection[0]);
				e0 = intersection[0];
				m_edgeBoundary[e0] = i;
				nb_boundary->set(nStart.id(), i);
				nb_boundary->set(nEnd.id(), i);
				nStart = e0.getOppositeNode(nStart);
			}
		}
		auto it = std::find(FrontEdges.begin(), FrontEdges.end(), e0);
		if (it != FrontEdges.end()) {
			// Si l'élément est trouvé, le supprimer
			FrontEdges.erase(it);
		}
		std::get<1>(m_nodeBoundary.at(nEnd))=e0;
		i++;
	}


	std::cout<<"Nombre de bords : "<<m_edgeBoundary.size()<<std::endl;
	std::cout<<"============== End Initialisation ==============\n"<<std::endl;
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("/home/pagea/Documents/Travail/data/toto_link.vtk");

}

/*----------------------------------------------------------------------------*/

Quadfront::STATUS
Quadfront::execute()
{
	return Quadfront::SUCCESS;
}

/*----------------------------------------------------------------------------*/

void Quadfront::set_edgeBoundary(Edge& Aedge, int i ){
	m_edgeBoundary[Aedge] = i;
}

/*----------------------------------------------------------------------------*/

int Quadfront::get_edgeBoundary(Edge& Aedge){
	return m_edgeBoundary[Aedge];
}

/*----------------------------------------------------------------------------*/

bool Quadfront::find_edgeBoundary(Edge& Aedge){
   return m_edgeBoundary.find(Aedge) != m_edgeBoundary.end();
}
/*----------------------------------------------------------------------------*/

void Quadfront::set_nodeBoundary(Node& Anode, int& i, Edge& Aedge0, Edge& Aedge1){
   m_nodeBoundary[Anode] = std::make_tuple(i, Aedge0, Aedge1);
}

/*----------------------------------------------------------------------------*/

std::tuple<int, Edge, Edge> Quadfront::get_nodeBoundary(Node& Anode){
	return m_nodeBoundary[Anode];
}

/*----------------------------------------------------------------------------*/

bool Quadfront::find_nodeBoundary(Node& Anode){
	return m_nodeBoundary.find(Anode) != m_nodeBoundary.end();
}

/*----------------------------------------------------------------------------*/

std::map<Node, std::pair<Edge,Edge>> Quadfront::get_edgeInFace(Face Aface){
	std::map<Node, std::pair<Edge,Edge>> result;
	std::vector<Edge> vEdge = Aface.get<Edge>();
	for (auto& n : Aface.get<Node>()){
		std::vector<Edge> vInter;
		std::vector<Edge> edgeNode = n.get<Edge>();
		std::sort(vEdge.begin(), vEdge.end());
		std::sort(edgeNode.begin(), edgeNode.end());
		std::set_intersection(vEdge.begin(), vEdge.end(), edgeNode.begin(), edgeNode.end(), std::back_inserter(vInter));
		result[n] = std::make_pair(vInter[0], vInter[1]);
	}
	return result;
}

/*----------------------------------------------------------------------------*/

Node Quadfront::get_opposit2Edge(Edge& Aedge, Face& Aface){
	std::vector<Node> vNode = Aface.get<Node>();
	for (auto& n : vNode) {
		if (n != Aedge.get<Node>()[0] and n != Aedge.get<Node>()[1]) {
			return n;
		}
	}
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::get_opposit2Node(Node& Anode, Face& Aface){
	std::vector<Edge> vEdge = Aface.get<Edge>();
	for (auto& e : vEdge) {
		if (Anode != e.template get<Node>()[0] and Anode != e.template get<Node>()[1]) {
			return e;
		}
	}
}

/*----------------------------------------------------------------------------*/

math::Vector3d Quadfront::node2Vector(Edge &Aedge, Node& Anode){
	return {Aedge.getOppositeNode(Anode).X() - Anode.X(),
	        Aedge.getOppositeNode(Anode).Y() - Anode.Y(),
	        Aedge.getOppositeNode(Anode).Z() - Anode.Z()};
}

/*----------------------------------------------------------------------------*/

double Quadfront::angle(Node Anode){
	// Calcul la somme des angles des triangles ayant Anode en commun
	std::vector<Face> vFace = Anode.get<Face>();
	double m_sommeAngle;
	for (auto& face : vFace){
		if (face.nbEdges() == 3){
			std::map<Node, std::pair<Edge,Edge>> node2edge = get_edgeInFace(face);
			m_sommeAngle += node2Vector(std::get<0>(node2edge[Anode]), Anode).angle(node2Vector(std::get<1>(node2edge[Anode]), Anode));
		}
	}
	return m_sommeAngle;
}

/*----------------------------------------------------------------------------*/

math::Vector3d Quadfront::get_bissectrice(math::Vector3d Avec1, math::Vector3d Avec2){
	if (Avec1.dot(Avec2)!=-1){
		return (Avec1 + Avec2).normalize();
	}
	else{
		return {Avec1[1],-Avec1[0],Avec1[2]};
	}
}

/*----------------------------------------------------------------------------*/

void Quadfront::initialise_nodeBoundary(){
	for (auto& n : m_nodeBoundary){
		if (angle(n.first)<=3*M_PI/4){
			std::get<0>(m_nodeBoundary[n.first]) = 1;
		}
	}
}

/*----------------------------------------------------------------------------*/

std::tuple<double, double, double> Quadfront::intersectionVec2Edge(math::Vector3d& vec_bissectrice, Node& Anode, Edge& Aedge){
	// a*Y+b*X+c = 0
	//Anode origine du vecteur vec
	math::Vector3d vecEdge = node2Vector(Aedge, Aedge.get<Node>()[0]);
	if (vec_bissectrice.X()==0 and Anode.X()==0){
		double vectorCoefc = -(vec_bissectrice.Y()*Anode.Y()+vec_bissectrice.Z()*Anode.Z());

		double edgeCoefc = -(vecEdge.Y()*Aedge.get<Node>()[0].Y()+
		                     vecEdge.Z()*Aedge.get<Node>()[0].Z());

		double v1 = -(vectorCoefc+vec_bissectrice.Y()*edgeCoefc/vecEdge.Y())/(vec_bissectrice.Z()-vec_bissectrice.Y()*vecEdge.Z()/vecEdge.Y());
		double v2 = (edgeCoefc-vecEdge.Z()*v1)/vecEdge.Y();
		return std::make_tuple(0, v1, v2);
	}
	else if (vec_bissectrice.Y()==0 and Anode.Y()==0){
		double vectorCoefc = -(vec_bissectrice.X()*Anode.X()+vec_bissectrice.Z()*Anode.Z());
		double edgeCoefc = -(vecEdge.X()*Aedge.get<Node>()[0].X()+
		                     vecEdge.Z()*Aedge.get<Node>()[0].Z());

		double v1 = -(vectorCoefc+vec_bissectrice.X()*edgeCoefc/vecEdge.X())/(vec_bissectrice.Z()-vec_bissectrice.X()*vecEdge.Z()/vecEdge.X());
		double v2 = (edgeCoefc-vecEdge.Z()*v1)/vecEdge.X();
		return std::make_tuple(v1, 0, v2);
	}
	else {
		double vectorCoefc = -(vec_bissectrice.X()*Anode.X()+vec_bissectrice.Y()*Anode.Y());
		double edgeCoefc = -(vecEdge.X()*Aedge.get<Node>()[0].X()+
		                     vecEdge.Y()*Aedge.get<Node>()[0].Y());

		double v1 = -(vectorCoefc+vec_bissectrice.X()*edgeCoefc/vecEdge.X())/(vec_bissectrice.Y()-vec_bissectrice.X()*vecEdge.Y()/vecEdge.X());
		double v2 = (edgeCoefc-vecEdge.Y()*v1)/vecEdge.X();
		return std::make_tuple(v1, v2, 0);
	}
}
/*----------------------------------------------------------------------------*/

void Quadfront::removeFace(Face& Aface){
	for (auto& n : Aface.get<Node>()){
		n.remove<Face>(Aface);
	}
	for (auto& e : Aface.get<Edge>()){
		e.remove<Face>(Aface);
	}
	m_mesh->deleteFace(Aface.id());
}

/*----------------------------------------------------------------------------*/

void Quadfront::removeEdge(Edge& Aedge){
	std::vector<Node> vNode = Aedge.get<Node>();
	vNode[0].remove<Edge>(Aedge);
	vNode[1].remove<Edge>(Aedge);
	for (auto& f : Aedge.get<Face>()){
		f.remove(Aedge);
	}
	m_mesh->deleteEdge(Aedge.id());
}

/*----------------------------------------------------------------------------*/


Edge Quadfront::swap(Edge& Aedge_inter,
                Node& Anode,
                Node& Anode_opposit,
                std::map<Node, std::pair<Edge, Edge>>& node2edge_opposit,
                Edge& Aedge_min1,
                Edge& Aedge_min2){

	std::vector<Face> face_inter = Aedge_inter.get<Face>();

	Node node1 = Aedge_min1.getOppositeNode(Anode);
	Node node2 = Aedge_min2.getOppositeNode(Anode);

	// Supprimer les faces
	std::cout<<"\n======Face supprimée 1 ======\n"<<std::endl;
	std::cout<<face_inter[0]<<std::endl;
	std::cout<<"\t==Edge dans la Face =="<<std::endl;
	for(auto& i2 : face_inter[0].get<Edge>()){
		std::cout<<"\t \tEdge "<<i2.id()<<" : "<<i2.get<Node>()[0]<<" ; "<<i2.get<Node>()[1]<<std::endl;
	}

	std::cout<<"\n======Face supprimée 2 ======\n"<<std::endl;
	std::cout<<face_inter[1]<<std::endl;
	std::cout<<"\t==Edge dans la Face =="<<std::endl;
	for(auto& x2 : face_inter[1].get<Edge>()){
		std::cout<<"\t \tEdge "<<x2.id()<<" : "<<x2.get<Node>()[0]<<" ; "<<x2.get<Node>()[1]<<std::endl;
	}

	removeFace(face_inter[0]);
	removeFace(face_inter[1]);


	std::cout<<"exemple 0 " <<face_inter[0]<<std::endl;
	std::cout<<"exemple 1 " <<face_inter[1]<<std::endl;

	// Création du de la diagonale
	Edge sideEdge = m_mesh->newEdge(Anode, Anode_opposit);


	Anode.add<Edge>(sideEdge);
	Anode_opposit.add<Edge>(sideEdge);

	// Suppression de l'arete opposée à Anode
	std::cout<<"\nEdge remove : "<<Aedge_inter.id()<<std::endl;
	removeEdge(Aedge_inter);

	// Création des deux nouveaux triangle
	Face newFace1 = m_mesh->newTriangle(Anode, Anode_opposit, node1);
	Face newFace2 = m_mesh->newTriangle(Anode, Anode_opposit, node2);

	// mettre à jour les connections des noeuds avec les faces
	Anode.add<Face>(newFace1);
	Anode.add<Face>(newFace2);
	Anode_opposit.add<Face>(newFace1);
	Anode_opposit.add<Face>(newFace2);
	node1.add<Face>(newFace1);
	node2.add<Face>(newFace2);

	// mettre à jour la connection des faces avec les edges
	newFace1.add<Edge>(sideEdge);
	newFace1.add<Edge>(Aedge_min1);
	newFace2.add<Edge>(sideEdge);
	newFace2.add<Edge>(Aedge_min2);

	sideEdge.add<Face>(newFace1);
	Aedge_min1.add<Face>(newFace1);
	sideEdge.add<Face>(newFace2);
	Aedge_min2.add<Face>(newFace2);

	//On cherche le coté opposé à edge_min
	if (std::get<0>(node2edge_opposit[node1]) != Aedge_inter) {
		newFace1.add<Edge>(std::get<0>(node2edge_opposit[node1]));
		std::get<0>(node2edge_opposit[node1]).add<Face>(newFace1);
	}
	else{
		newFace1.add<Edge>(std::get<1>(node2edge_opposit[node1]));
		std::get<1>(node2edge_opposit[node1]).add<Face>(newFace1);
	}
	if (std::get<0>(node2edge_opposit[node2]) != Aedge_inter) {
		newFace2.add<Edge>(std::get<0>(node2edge_opposit[node2]));
		std::get<0>(node2edge_opposit[node2]).add<Face>(newFace2);
	}
	else{
		newFace2.add<Edge>(std::get<1>(node2edge_opposit[node2]));
		std::get<1>(node2edge_opposit[node2]).add<Face>(newFace2);
	}


	std::cout<<"\n==== Nouvelle Face 1 ====\n"<<std::endl;
	std::cout<<newFace1<<std::endl;
	std::cout<<"\t==Edge dans la Face =="<<std::endl;
	for(auto& e1 : newFace1.get<Edge>()){
		std::cout<<"\t \tEdge "<<e1.id()<<" : "<<e1.get<Node>()[0]<<" ; "<<e1.get<Node>()[1]<<std::endl;
		std::cout<<"\t \tEdge size "<<e1.get<Face>().size()<<std::endl;
	}

	std::cout<<"\n==== Nouvelle Face 2 ====\n"<<std::endl;
	std::cout<<newFace2<<std::endl;
	std::cout<<"\t==Edge dans la Face =="<<std::endl;
	for(auto& e2 : newFace2.get<Edge>()){
		std::cout<<"\t \tEdge "<<e2.id()<<" : "<<e2.get<Node>()[0]<<" ; "<<e2.get<Node>()[1]<<std::endl;
		std::cout<<"\t \tEdge size "<<e2.get<Face>().size()<<std::endl;
	}

	std::cout<<"\n"<<std::endl;

return sideEdge;
}

/*----------------------------------------------------------------------------*/

Edge Quadfront::split(Node& AnewNode, Edge& Aedge_inter,
                 Node& Anode,
                 Node& Anode_opposit,
                 std::map<Node, std::pair<Edge, Edge>>& node2edge_opposit,
                 Edge& Aedge_min1,
                 Edge& Aedge_min2){

	Node node1 = Aedge_min1.getOppositeNode(Anode);
	Node node2 = Aedge_min2.getOppositeNode(Anode);


	std::vector<Face> face_inter = Aedge_inter.get<Face>();

	// Création du de la diagonale et actualisation des connections avec les noeuds
	Edge sideEdge = m_mesh->newEdge(Anode, AnewNode);
	Edge sideEdge_opposit = m_mesh->newEdge(AnewNode, Anode_opposit);

	Anode.add<Edge>(sideEdge);
	AnewNode.add<Edge>(sideEdge);
	AnewNode.add<Edge>(sideEdge_opposit);
	Anode_opposit.add<Edge>(sideEdge_opposit);

	// CHECKER LE SENS DES EDGES
	Edge splitEdge1 = m_mesh->newEdge(node1, AnewNode);
	Edge splitEdge2 = m_mesh->newEdge(AnewNode, node2);

	node1.add<Edge>(splitEdge1);
	node2.add<Edge>(splitEdge2);
	AnewNode.add<Edge>(splitEdge1);
	AnewNode.add<Edge>(splitEdge2);

	// Création des deux nouveaux triangle
	Face newFace1 = m_mesh->newTriangle(Anode, AnewNode, node1);
	Face newFace2 = m_mesh->newTriangle(Anode, AnewNode, node2);
	Face newFace_opposit1 = m_mesh->newTriangle(Anode_opposit, AnewNode, node1);
	Face newFace_opposit2 = m_mesh->newTriangle(Anode_opposit, AnewNode, node2);


	// mettre à jour les connections des noeuds avec les faces
	Anode.add<Face>(newFace1);
	Anode.add<Face>(newFace2);
	AnewNode.add<Face>(newFace1);
	AnewNode.add<Face>(newFace2);
	node1.add<Face>(newFace1);
	node2.add<Face>(newFace2);

	Anode_opposit.add<Face>(newFace_opposit1);
	Anode_opposit.add<Face>(newFace_opposit2);
	AnewNode.add<Face>(newFace_opposit1);
	AnewNode.add<Face>(newFace_opposit2);
	node1.add<Face>(newFace_opposit1);
	node2.add<Face>(newFace_opposit2);

	// mettre à jour la connection des faces avec les edges
	//CHECKER
	newFace1.add<Edge>(sideEdge);
	newFace2.add<Edge>(sideEdge);
	newFace_opposit1.add<Edge>(sideEdge_opposit);
	newFace_opposit2.add<Edge>(sideEdge_opposit);
	newFace1.add<Edge>(splitEdge1);
	newFace2.add<Edge>(splitEdge2);
	newFace_opposit1.add<Edge>(splitEdge1);
	newFace_opposit2.add<Edge>(splitEdge2);
	newFace1.add<Edge>(Aedge_min1);
	newFace2.add<Edge>(Aedge_min2);

	sideEdge.add<Face>(newFace1);
	sideEdge.add<Face>(newFace2);
	sideEdge_opposit.add<Face>(newFace_opposit1);
	sideEdge_opposit.add<Face>(newFace_opposit2);
	splitEdge1.add<Face>(newFace_opposit1);
	splitEdge1.add<Face>(newFace1);
	splitEdge2.add<Face>(newFace_opposit2);
	splitEdge2.add<Face>(newFace2);
	Aedge_min1.add<Face>(newFace1);
	Aedge_min2.add<Face>(newFace2);

	//On cherche le coté opposé à edge_min
	if (std::get<0>(node2edge_opposit[node1]) != Aedge_inter) {
		newFace_opposit1.add<Edge>(std::get<0>(node2edge_opposit[node1]));
		std::get<0>(node2edge_opposit[node1]).add<Face>(newFace_opposit1);
	}
	else{
		newFace_opposit1.add<Edge>(std::get<1>(node2edge_opposit[node1]));
		std::get<1>(node2edge_opposit[node1]).add<Face>(newFace_opposit1);
	}
	if (std::get<0>(node2edge_opposit[node2]) != Aedge_inter) {
		newFace_opposit2.add<Edge>(std::get<0>(node2edge_opposit[node1]));
		std::get<0>(node2edge_opposit[node1]).add<Face>(newFace_opposit2);
	}
	else{
		newFace_opposit2.add<Edge>(std::get<1>(node2edge_opposit[node1]));
		std::get<1>(node2edge_opposit[node1]).add<Face>(newFace_opposit2);
	}


	// Suppression de l'arete opposée à Anode
	removeEdge(Aedge_inter);

	// Supprimer les faces
	removeFace(face_inter[0]);
	removeFace(face_inter[1]);

	return sideEdge;
}

/*----------------------------------------------------------------------------*/

// checker au préalable le statut du noeud => dans le cas ou c'est 0
// < pi/6 pour que ce soit valable
Edge Quadfront::sideEdge(Node& Anode, Edge& Aedge){

	Edge sideEdge;

	if (std::get<0>(m_nodeBoundary[Anode])==1){
		if (std::get<1>(m_nodeBoundary[Anode])==Aedge){
			return std::get<2>(m_nodeBoundary[Anode]);
		}
		else{
			return std::get<1>(m_nodeBoundary[Anode]);
		}
	}

	Edge a = std::get<1>(m_nodeBoundary[Anode]);
	Edge b = std::get<2>(m_nodeBoundary[Anode]);
	math::Vector3d vE0 = node2Vector(a, Anode).normalize();
	math::Vector3d vE1 = node2Vector(b, Anode).normalize();


	// Bissectrice des arêtes adjacentes au noeud
	math::Vector3d m_bissectrice = get_bissectrice(vE0, vE1).normalize();
	double dot_max1 = 0;
	double dot_max2 = 0;
	Edge edge_min1;
	Edge edge_min2;

	for (Edge& e : Anode.get<Edge>()) {
		math::Vector3d vec = node2Vector(e, Anode).normalize();
		if (m_bissectrice.dot(vec) > dot_max1) {
			dot_max1 = m_bissectrice.dot(vec);
			edge_min1 = e;
		}
		else if (m_bissectrice.dot(vec) > dot_max2) {
			dot_max2 = m_bissectrice.dot(vec);
			edge_min2 = e;
		}
	}
	std::cout<<"Edge min 1 : "<< edge_min1.get<Node>()[0]<<"\t"<<edge_min1.get<Node>()[1]<<std::endl;
	std::cout<<"Edge min 2 : "<< edge_min2.get<Node>()[0]<<"\t"<<edge_min2.get<Node>()[1]<<std::endl;

	//CAS IDEAL
	if (acos(dot_max1) < M_PI / 24) {
		sideEdge = edge_min1;
		std::cout<<"l'angle "<<edge_min1.id()<<" est plus petit que Pi/6 "<<std::endl;
	}
	else {
		Node node1 = edge_min1.getOppositeNode(Anode);
		Node node2 = edge_min2.getOppositeNode(Anode);

		std::vector<TCellID> vFace = m_mesh->getCommonFaces(node1, node2);

		std::map<Node, std::pair<Edge, Edge>> node2edge;
		std::map<Node, std::pair<Edge, Edge>> node2edge_opposit;
		Node opposit_node;
		Face face;
		Face opposit_face;

		if (get_edgeInFace(m_mesh->get<Face>(vFace[0])).find(Anode) != get_edgeInFace(m_mesh->get<Face>(vFace[0])).end()){
			face = m_mesh->get<Face>(vFace[0]);
			opposit_face = m_mesh->get<Face>(vFace[1]);
			node2edge = get_edgeInFace(m_mesh->get<Face>(vFace[0]));
			node2edge_opposit = get_edgeInFace(m_mesh->get<Face>(vFace[1]));
		}
		else{
			face = m_mesh->get<Face>(vFace[1]);
			opposit_face = m_mesh->get<Face>(vFace[0]);
			node2edge = get_edgeInFace(m_mesh->get<Face>(vFace[1]));
			node2edge_opposit = get_edgeInFace(m_mesh->get<Face>(vFace[0]));
		}

		for (auto& n : node2edge_opposit){
			if (n.first != node1 and n.first != node2){
				opposit_node = n.first;
			}
		}

		math::Vector3d vec = {opposit_node.X() - Anode.X(), opposit_node.Y() - Anode.Y(), opposit_node.Z() - Anode.Z()};

		std::vector<Edge> vFace_ = face.get<Edge>();
		std::vector<Edge> vFace_opposit= opposit_face.get<Edge>();
		std::vector<Edge> vFace_inter;

		std::set_intersection(vFace_.begin(), vFace_.end(), vFace_opposit.begin(), vFace_opposit.end(), std::back_inserter(vFace_inter));
		Edge edge_inter = vFace_inter[0];


		// SWAP CASE
		if (m_bissectrice.angle(vec) < M_PI / 6 and vec.norm() < sqrt(3) * (edge_min1.length() + edge_min2.length()) / 2) {
			std::cout<<"==== SWAP ===="<<std::endl;
			sideEdge = swap(edge_inter, Anode, opposit_node, node2edge_opposit, edge_min1, edge_min2);
		}

		// SPLIT CASE
		else {
			std::tuple<double, double, double> vNode = intersectionVec2Edge(m_bissectrice, Anode, edge_inter);
			Node newNode = m_mesh->newNode(std::get<0>(vNode), std::get<1>(vNode), std::get<2>(vNode));

			std::cout<<"==== SPLIT ===="<<std::endl;
			sideEdge = split(newNode, edge_inter, Anode, opposit_node, node2edge_opposit,edge_min1, edge_min2);

		}
	}
	return sideEdge;
}
/*----------------------------------------------------------------------------*/

std::vector<Face> Quadfront::get_superposed(Node& Anode, Edge& Aedge, Edge& sideEdge0, Edge& sideEdge1, Edge& upEdge){
	std::vector<Face> vFace_in;


}


std::vector<Edge>  Quadfront::interestionS(Node& AnodeRecovery0, Node& AnodeRecovery1){

	Variable<int>* sequent = m_mesh->newVariable<int, GMDS_EDGE>("sequent");

	math::Vector3d vecAB({AnodeRecovery1.X()-AnodeRecovery0.X(),
	                      AnodeRecovery1.Y()-AnodeRecovery0.Y(),
	                      AnodeRecovery1.Z()-AnodeRecovery0.Z()});


	math::Vector3d vecNormaleAB({vecAB[1],-vecAB[0], vecAB[2]});


	// CHECKER POUR LA 3D
	std::vector<Face> vFace0 = AnodeRecovery0.get<Face>();
	std::vector<Edge> vEdgeIntersect;
	std::map<Node, std::pair<Edge, Edge>> node2edgeInFace;
	Edge e_opposit;
	Face face;


	//determine la meilleure face (celle dont une arete se situe en dessous et l'autre au dessus)
	for (auto& f : vFace0){

		node2edgeInFace = get_edgeInFace(f);

		math::Vector3d vec0 = node2Vector(std::get<0>(node2edgeInFace[AnodeRecovery0]), AnodeRecovery0);
		math::Vector3d vecNormal0 = {vec0[1], -vec0[0], vec0[2]};
		math::Vector3d vec1 = node2Vector(std::get<1>(node2edgeInFace[AnodeRecovery0]), AnodeRecovery0);
		math::Vector3d vecNormal1 = {vec1[1], -vec1[0], vec0[2]};




		// Dans le cas ou une arête est alignée avec le segment AB
		if (vecAB.dot(vec0)==1){
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).X() += 0.05 * vecNormaleAB.X();
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Y() += 0.05 * vecNormaleAB.Y();
			std::get<0>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Z() += 0.05 * vecNormaleAB.Z();
			vec0 = vec0 + 0.05 * vecNormaleAB;
			vecNormal0 = {vec0[1], -vec0[0], vec0[2]};
		}
		if (vecAB.dot(vec1)==1){
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).X() += 0.05 * vecNormaleAB.X();
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Y() += 0.05 * vecNormaleAB.Y();
			std::get<1>(node2edgeInFace[AnodeRecovery0]).getOppositeNode(AnodeRecovery0).Z() += 0.05 * vecNormaleAB.Z();
			vec1 = vec1 + 0.05 * vecNormaleAB;
			vecNormal1 = {vec1[1], -vec1[0], vec1[2]};

		}
		if ((vecAB.dot(vecNormal0)>0 and vecAB.dot(vecNormal1)<0) or (vecAB.dot(vecNormal0)<0 and vecAB.dot(vecNormal1)>0)){
			if (vec0.dot(vecAB)>0 and vec1.dot(vecAB)>0){
				e_opposit = get_opposit2Node(AnodeRecovery0, f);
				face = f;
			}
		}
	}

	int i =0;

	while (m_edgeBoundary.find(e_opposit) == m_edgeBoundary.end()){
		i++;
		vEdgeIntersect.push_back(e_opposit);
		sequent->set(e_opposit.id(), i);
		std::cout<<"e_opposit "<<e_opposit.id()<<" : "<<e_opposit.get<Node>()[0]<<" ; "<<e_opposit.get<Node>()[1]<<std::endl;

		if(e_opposit.get<Face>()[0] == face){
			face = e_opposit.get<Face>()[1];
		}
		else{
			face = e_opposit.get<Face>()[0];
		}


		node2edgeInFace = get_edgeInFace(face);

		if (node2edgeInFace.find(AnodeRecovery1) == node2edgeInFace.end()){
			Node n_opposit = get_opposit2Edge(e_opposit, face);

			std::cout<<"Node opposite "<<n_opposit<<std::endl;

			math::Vector3d vecAC = {n_opposit.X()-AnodeRecovery0.X(),
			                        n_opposit.Y()-AnodeRecovery0.Y(),
			                        n_opposit.Z()-AnodeRecovery0.Z()};


			double test = 1 - vecAB.dot(vecAC);

			if (test == 0.){
				if (m_nodeBoundary.find(n_opposit) == m_nodeBoundary.end()){
					std::cout<<"Node avant modif"<<n_opposit<<std::endl;
					n_opposit.X() += 0.05 * vecNormaleAB.X();
					n_opposit.Y() += 0.05 * vecNormaleAB.Y();
					n_opposit.Z() += 0.05 * vecNormaleAB.Z();
					vecAC = vecAC + 0.05*vecNormaleAB;
					std::cout<<"Node après modif "<<n_opposit<<std::endl;
				}
			}
			math::Vector3d vecNormaleAC = {vecAC[1], -vecAC[0], vecAC[2]};

			Node node0 = std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit);

			math::Vector3d vec0 = {std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).X()-AnodeRecovery0.X(),
			                       std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Y()-AnodeRecovery0.Y(),
			                       std::get<0>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Z()-AnodeRecovery0.Z()};


			math::Vector3d vecNormale0 = {vec0[1], -vec0[0], vec0[2]};

			math::Vector3d vec1 = {std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).X()-AnodeRecovery0.X(),
			                       std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Y()-AnodeRecovery0.Y(),
			                       std::get<1>(node2edgeInFace[n_opposit]).getOppositeNode(n_opposit).Z()-AnodeRecovery0.Z()};


			math::Vector3d vecNormale1 = {vec0[1], -vec0[0], vec0[2]};



			if (vecAB.dot(vecNormaleAC)>0){
				if (vecAB.dot(vecNormale0)>0){
					e_opposit = std::get<1>(node2edgeInFace[n_opposit]);
				}
				else{
					e_opposit = std::get<0>(node2edgeInFace[n_opposit]);
				}
			}
			else{
				if (vecAB.dot(vecNormale0)>0){
					e_opposit = std::get<0>(node2edgeInFace[n_opposit]);
				}
				else{
					e_opposit = std::get<1>(node2edgeInFace[n_opposit]);
				}
			}
		}
		else{
			break;
		}
	}
	return vEdgeIntersect;
}

/*----------------------------------------------------------------------------*/

void Quadfront::testInside(Face& faceRef, Face& faceInside, std::vector<Face>& listInside){
	std::vector<Edge> edgeFaceRef = faceRef.get<Edge>();

	if (std::find(listInside.begin(), listInside.end(), faceInside) != listInside.end() or faceInside.nbEdges() == 4){
		return;
	}
	bool test = true;
	double Xmax = std::max(std::max(faceRef.get<Node>()[0].X(),faceRef.get<Node>()[1].X()), std::max(faceRef.get<Node>()[2].X(), faceRef.get<Node>()[3].X()));
	double Xmin = std::min(std::min(faceRef.get<Node>()[0].X(),faceRef.get<Node>()[1].X()), std::min(faceRef.get<Node>()[2].X(), faceRef.get<Node>()[3].X()));
	double Ymax = std::max(std::max(faceRef.get<Node>()[0].Y(),faceRef.get<Node>()[1].Y()), std::max(faceRef.get<Node>()[2].Y(), faceRef.get<Node>()[3].Y()));
	double Ymin = std::min(std::min(faceRef.get<Node>()[0].Y(),faceRef.get<Node>()[1].Y()), std::min(faceRef.get<Node>()[2].Y(), faceRef.get<Node>()[3].Y()));
	double Zmax = std::max(std::max(faceRef.get<Node>()[0].Z(),faceRef.get<Node>()[1].Z()), std::max(faceRef.get<Node>()[2].Z(), faceRef.get<Node>()[3].Z()));
	double Zmin = std::min(std::min(faceRef.get<Node>()[0].Z(),faceRef.get<Node>()[1].Z()), std::min(faceRef.get<Node>()[2].Z(), faceRef.get<Node>()[3].Z()));
	for (auto& n : faceInside.get<Node>()){
		if ((n.X()<Xmin or n.X()>Xmax) or (n.Y()<Ymin or n.Y()>Ymax) or (n.Z()<Zmin or n.Z()>Zmax)){
			test = false;
		}
	}
	if (test){
		listInside.push_back(faceInside);
		for (auto& e : faceInside.get<Edge>()){
			if (std::find(edgeFaceRef.begin(), edgeFaceRef.end(), e) == edgeFaceRef.end()){
				if (e.get<Face>()[0] != faceInside){
					testInside(faceRef, e.get<Face>()[0], listInside);
				}
				else{
					testInside(faceRef, e.get<Face>()[1], listInside);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/

std::vector<Face> Quadfront::createQuad(Node& Anode, Edge& Aedge, Edge& sideEdge0, Edge& sideEdge1, Edge& recoveryEdge){
	Face face0;
	if (sideEdge0.get<Face>()[0].nbEdges() == 3){
		face0 = sideEdge0.get<Face>()[0];
		std::cout<<" face0 : "<<face0<<std::endl;
	}
	else{
		face0 = sideEdge0.get<Face>()[1];
		std::cout<<" face0 : "<<face0<<std::endl;
	}
	std::vector<Face> listInside;

	Face quad = m_mesh->newQuad(Aedge.get<Node>()[0],
	   sideEdge0.getOppositeNode(Aedge.get<Node>()[0]),
	                            Aedge.get<Node>()[1], sideEdge1.getOppositeNode(Aedge.get<Node>()[1]));

	quad.add<Edge>(Aedge);
	quad.add<Edge>(sideEdge0);
	quad.add<Edge>(sideEdge1);
	quad.add<Edge>(recoveryEdge);



	testInside(quad, face0, listInside);

	Aedge.add<Face>(quad);
	sideEdge0.add<Face>(quad);
	sideEdge1.add<Face>(quad);
	recoveryEdge.add<Face>(quad);

	for (auto& f : listInside){
		removeFace(f);
	}

	for (auto& f : listInside){
		std::cout<<f<<std::endl;
	}

	return listInside;
}

/*----------------------------------------------------------------------------*/

void  Quadfront::edgeRecovery(Edge Aedge){
	Node node0 =Aedge.get<Node>()[0];
	Edge sideEdge0;
	Node nodeRecovery0;
	Edge edgeRecovery;

	std::cout<<"\n======= SideEdge0 ========\n"<<std::endl;
	if (std::get<0>(m_nodeBoundary[node0])==0){
		sideEdge0 = sideEdge(node0, Aedge);
		nodeRecovery0 = sideEdge0.getOppositeNode(node0);
	}
	else{
		if (std::get<1>(m_nodeBoundary[node0])==Aedge){
			sideEdge0 = std::get<2>(m_nodeBoundary[node0]);
			nodeRecovery0 = sideEdge0.getOppositeNode(node0);
		}
		else{
			sideEdge0 = std::get<1>(m_nodeBoundary[node0]);
			nodeRecovery0 = sideEdge0.getOppositeNode(node0);
		}
	}
	std::cout<<"sideEdge0 "<<sideEdge0.id()<<" : "<<sideEdge0.get<Node>()[0]<<" ; "<<sideEdge0.get<Node>()[1]<<std::endl;
	std::cout<<"nodeRecovery0 : "<<nodeRecovery0<<std::endl;

	std::cout<<"\n======= SideEdge1 ========\n"<<std::endl;
	Node node1 =Aedge.get<Node>()[1];
	Edge sideEdge1;
	Node nodeRecovery1;
	if (std::get<0>(m_nodeBoundary[node1])==0){
		sideEdge1 = sideEdge(node1, Aedge);
		nodeRecovery1 = sideEdge1.getOppositeNode(node1);
	}
	else{
		if (std::get<1>(m_nodeBoundary[node1])==Aedge){
			sideEdge1 = std::get<2>(m_nodeBoundary[node1]);
			nodeRecovery1 = sideEdge1.getOppositeNode(node1);
		}
		else{
			sideEdge1 = std::get<1>(m_nodeBoundary[node1]);
			nodeRecovery1 = sideEdge1.getOppositeNode(node1);
		}
	}
	std::cout<<"sideEdge1 "<<sideEdge1.id()<<" : "<<sideEdge1.get<Node>()[0]<<" ; "<<sideEdge1.get<Node>()[1]<<std::endl;
	std::cout<<"nodeRecovery1 : "<<nodeRecovery1<<std::endl;



	// CAS PARTICULIER nodeRecovery0 = nodeRecovery1
	//if (nodeRecovery0 == nodeRecovery1){
	//}
	//  A FAIRE


	math::Vector3d vecAB({nodeRecovery1.X()-nodeRecovery0.X(),
	                      nodeRecovery1.Y()-nodeRecovery0.Y(),
	                      nodeRecovery1.Z()-nodeRecovery0.Z()});

	std::vector<Face> Face0 = sideEdge0.get<Face>();
	std::vector<Face> Face1 = sideEdge1.get<Face>();
	std::vector<Face> interFace;

	std::cout<<"\n ==== UpSide ====\n "<<std::endl;

	std::vector<Edge> listIntersect = interestionS(nodeRecovery0, nodeRecovery1);

	for (auto e = listIntersect.begin(); e!=listIntersect.end();){
		std::vector<Node> vNode = (*e).get<Node>();
		std::vector<Face> vFace = (*e).get<Face>();
		Node n_opposit0 = get_opposit2Edge(*e, vFace[0]);
		std::cout<<"n_opposit0 : "<<n_opposit0<<std::endl;
		std::map<Node, std::pair<Edge, Edge>> node2edgeInFace0 = get_edgeInFace(vFace[0]);
		Node n_opposit1 = get_opposit2Edge(*e, vFace[1]);
		std::cout<<"n_opposit1 : "<<n_opposit1<<std::endl;
		std::map<Node, std::pair<Edge, Edge>> node2edgeInFace1 = get_edgeInFace(vFace[1]);
		math::Triangle triangle0(n_opposit0.point(), n_opposit1.point(), vNode[1].point());
		math::Triangle triangle1(n_opposit0.point(), n_opposit1.point(), vNode[0].point());


		if (triangle0.area()>0 and triangle1.area()>0){
			Edge diagEdge0;


			diagEdge0 = swap(*e, n_opposit0, n_opposit1, node2edgeInFace1, std::get<0>(node2edgeInFace0[n_opposit0]), std::get<1>(node2edgeInFace0[n_opposit0]));

			if ((diagEdge0.get<Node>()[0] == nodeRecovery0 and diagEdge0.get<Node>()[1] == nodeRecovery1) or
			    (diagEdge0.get<Node>()[1] == nodeRecovery0 and diagEdge0.get<Node>()[0] == nodeRecovery1)){
				edgeRecovery = diagEdge0;
			}

			std::cout<<"diagEdge0 : "<<diagEdge0.get<Node>()[0]<<diagEdge0.get<Node>()[1]<<std::endl;
			listIntersect.erase(e);

			std::tuple<double, double, double> vNode = intersectionVec2Edge(vecAB, nodeRecovery0, diagEdge0);
			math::Point point = (std::get<0>(vNode), std::get<1>(vNode), std::get<2>(vNode));
			if (diagEdge0.segment().isIn(point)){
				listIntersect.push_back(diagEdge0);
			}
		}
		else{
			listIntersect.erase(e);
			listIntersect.push_back(*e);
		}
	}


	std::cout<<"edgeRecovery : "<<edgeRecovery.get<Node>()[0]<<" ; "<<edgeRecovery.get<Node>()[1]<<std::endl;
	std::vector<Face> listInside = createQuad(node0, Aedge, sideEdge0, sideEdge1, edgeRecovery);
	std::cout<<" size : "<<listInside.size()<<std::endl;
}

/*----------------------------------------------------------------------------*/
}  // end namespace quadfront
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/