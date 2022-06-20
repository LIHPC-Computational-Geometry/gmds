//
// Created by Paul Bourmaud on 12/04/2022.
//

#include "gmds/paul/Tools.h"

using namespace gmds;
Tools::Tools(GridBuilderAround *AGrid)
	:g_grid(*AGrid){;}

int Tools::getValueActivateFace(const int faceIDChecked)
{
	int valueActivate;
	valueActivate = g_grid.activate->value(faceIDChecked);
	return valueActivate;

}

bool Tools::checkExistEdge(const int i1, const int i2, const int faceID)
{
	if(Tools::checkFollowIdNode(i1,i2,faceID) && Tools::checkCommonFace(i1,i2)){
		//std::cout<<"LES 2 NOEUDS ONT UNE ARETE !"<<std::endl;
		return true;
	}
	else{
		//std::cout<<"LES 2 NOEUDS N'ONT PAS UNE ARETE !"<<std::endl;
		return false;
	}

}

bool Tools::checkCommonFace(const int i1, const int i2)
{
	std::vector<Face> list_f= Tools::getListFacesOfNode(i1);
	for (auto f : list_f){
		std::vector<Node> list_n = Tools::getListNodesOfFace(f.id());
		for (auto n : list_n){
			if (n.id() == i2 ){
				//std::cout<<"LES 2 NOEUDS ONT UNE FACE COMMUNE"<<std::endl;
				return true;
			}
		}
	}
	//std::cout<<"LES 2 NOEUDS N'ONT PAS UNE FACE COMMUNE"<<std::endl;
	return false;
}

bool Tools::checkFollowIdNode(const int i1, const int i2,const int faceID)
{
	if (i1 == Tools::getIdNextNode(i2,faceID) || i1 == Tools::getIdPreviousNode(i2,faceID)){
		//std::cout<<"LES 2 NOEUDS SE SUIVENT !"<<std::endl;
		return true;
	}
	else{
		//std::cout<<"LES 2 NOEUDS NE SE SUIVENT PAS !"<<std::endl;
		return false;
	}
}

std::vector<Node> Tools::getListNodesOfFace(const int faceID)
{
	Face f = g_grid.m_mesh->get<Face>(faceID);
	std::vector<Node>list_nodes = f.get<Node>();
	/*
	for (auto n : list_nodes){
		std::cout<<"Face :"<< faceID <<" avec les noeuds :"<< n<<std::endl;
	}*/
	return list_nodes;
}

std::vector<Face> Tools::getListFacesOfNode(const int nodeID)
{
	Node n =g_grid.m_mesh->get<Node>(nodeID);
	std::vector<Face> list_faces = n.get<Face>();
	/*
	for (auto n : list_faces){
		std::cout<<"Le noeud " << nodeID << " avec les faces"<< n <<std::endl;
	}*/
	return list_faces;
}

std::vector<Face> Tools::getFacesCommon(const int i1, const int i2)
{
	std::vector<Face> list_Face_Common;
	std::vector<Face> list_f_i1 = Tools::getListFacesOfNode(i1);
	std::vector<Face> list_f_i2 = Tools::getListFacesOfNode(i2);
	for (auto f1 : list_f_i1){
		for (auto f2  : list_f_i2){
			if (f2 == f1 ){
				list_Face_Common.insert(list_Face_Common.end(),f1);
			}
		}
	}
	if(list_Face_Common.empty()){
		//std::cout<<"Erreur, noeuds pas de face en commun "<<std::endl;
		return {};
	}
	/*std::cout<<"List des faces en commun entre "<<i1<<" et "<<i2<<" :"<<std::endl;
	for (auto f : list_Face_Common){
		std::cout<< f <<std::endl;
	}*/
	return list_Face_Common;
}

int Tools::getIdOneCommonFace(const int i1, const int i2)
{
	std::vector<Face> list_f= Tools::getListFacesOfNode(i1);
	for (auto f : list_f){
		std::vector<Node> list_n = Tools::getListNodesOfFace(f.id());
		for (auto n : list_n){
			if (n.id() == i2 ){
				//std::cout<<"LES 2 NOEUDS ONT UNE FACE COMMUNE"<<std::endl;
				return f.id();
			}
		}
	}
	//std::cout<<"LES 2 NOEUDS N'ONT PAS UNE FACE COMMUNE"<<std::endl;
	return {};

}

int Tools::getIdPreviousNode(const int idNode,const int idFaceNode)
{
	std::vector<Node> listNodeFace = Tools::getListNodesOfFace(idFaceNode);

	for(int i = 0; i < listNodeFace.size();i++){
		//std::cout<<"L'element i : "<< i <<" du vecteur est :"<<listNodeFace[i].id()<<std::endl;
		//std::cout<<"Le dernier element du vec est "<<listNodeFace.back()<<std::endl;
		if (listNodeFace[i].id() == idNode){
			if(i == 0){
				//std::cout << "Le noeud precedent est " << listNodeFace.back().id() << std::endl;
				return listNodeFace.back().id();
			}
			else {
				//std::cout << "Le noeud precedent est " << listNodeFace[i - 1].id() << std::endl;
				return listNodeFace[i - 1].id();
			}
		}
	}
}

int Tools::getIdNextNode(const int idNode, const int idFaceNode)
{
	std::vector<Node> listNodeFace = Tools::getListNodesOfFace(idFaceNode);

	for(int i = 0; i < listNodeFace.size();i++){
		//std::cout<<"L'element i : "<< i <<" du vecteur est :"<<listNodeFace[i].id()<<std::endl;
		//std::cout<<"L'element idNode : "<< idNode <<std::endl;
		//std::cout<<"Le dernier element du vec est "<<listNodeFace.back()<<std::endl;
		if (listNodeFace[i].id() == idNode){
			if(i == listNodeFace.size()-1){
				//std::cout << "Le noeud suivant est " << listNodeFace.front().id() << std::endl;
				return listNodeFace.front().id();
			}
			else {
				//std::cout << "Le noeud suivant est " << listNodeFace[i + 1].id() << std::endl;
				return listNodeFace[i + 1].id();
			}
		}
	}
}

std::vector<std::vector<int>> Tools::getOtherNodes(const int i1, const int i2)
{
	std::vector<Face> listFace = Tools::getFacesCommon(i1,i2);
	std::vector<Node> listOtherNodes;
	std::vector<std::vector<int>> listPairNodes;
	std::vector<int> listElementAdd;

	for(auto f : listFace){
		std::vector<Node> nodesFace = Tools::getListNodesOfFace(f.id());
		for (auto n : nodesFace){
			if(n.id() != i1 && n.id() != i2){
				listOtherNodes.insert(listOtherNodes.end(),n);
			}
		}
	}

	for (auto n : listOtherNodes){
		for(auto n1 : listOtherNodes){
			if (n != n1 && Tools::checkCommonFace(n.id(),n1.id())){
				if (std::any_of(listElementAdd.begin(),listElementAdd.end(),[&](const int& elem) { return elem == n.id(); })
					 && std::any_of(listElementAdd.begin(),listElementAdd.end(),[&](const int& elem) { return elem == n1.id(); })){
					//std::cout<<"DEJA DANS LA LISTE :"<<n<< " & "<<n1<<std::endl;
				}
				else{
					std::vector<int> pairNodes;
					listElementAdd.push_back(n.id());
					listElementAdd.push_back(n1.id());
					pairNodes.insert(pairNodes.end(),n.id());
					pairNodes.insert(pairNodes.end(),n1.id());
					listPairNodes.push_back(pairNodes);
				}
			}
		}
	}
	/*for (auto n : listPairNodes){
		std::cout<<"UNE PAIRE de noeuds : "<<std::endl;
		for(auto i : n) {
			std::cout << "============"<< "Elements list pair : " << i << "===========" << std::endl;
		}
	}*/
	return listPairNodes;
}

std::vector<Face> Tools::getAllFacesChain(const int i1, const int i2)
{
	std::vector<Face> listAllFaces = {};
	std::vector<std::vector<int>> listAllPairNodes = Tools::getOtherNodes(i1,i2);
	std::vector<int> nodeCheck = {i1,i2};
	auto listAllPairToCheck = Tools::getAllNodesChain(i1,i2);
	auto listPairToCheck = listAllPairNodes;
	/*std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
	for (auto n : listAllPairToCheck){
		std::cout<<"une PAIR : "<<n.front()<<" "<<n.back()<<std::endl;
	}
	std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
	 */

	for (auto p : listAllPairToCheck){
		auto theCurrentFace = Tools::getFacesCommon(p.front(),p.back());
		for (auto f : theCurrentFace){
			if (std::find(listAllFaces.begin(), listAllFaces.end(), f)!=listAllFaces.end()){
				//std::cout<<"Other nodes deja check"<<std::endl;
			}
			else{
				listAllFaces.push_back(f);
			}
			//std::cout<<f<<std::endl;
		}
	}
	/*std::cout<<"LA LISTE DES FACES TRUE "<<std::endl;
	for(auto f : listAllFaces){
		std::cout<<f<<std::endl;

	}*/
	return listAllFaces;
}

std::vector<std::vector<int>> Tools::getAllNodesChain(const int i1, const int i2)
{
	std::vector<int> pairNodes={i1,i2};
	auto nodeCheck = pairNodes;
	std::vector<std::vector<int>> listAllPairNodesChain= {pairNodes};
	std::vector<std::vector<int>> listPairToCheck = Tools::getOtherNodes(i1,i2);

	while(!listPairToCheck.empty()) {
		auto listParcour = listPairToCheck;
		/*std::cout<<"List Parcours :"<<std::endl;
		for (auto n : listParcour){
			for(auto i : n){
				std::cout<<i<<std::endl;
			}
		}*/
		for (auto p : listParcour) {
			listAllPairNodesChain.push_back(p);
			nodeCheck.push_back(p.front());
			nodeCheck.push_back(p.back());
			auto newPairToAdd = Tools::getOtherNodes(p.front(), p.back());
			for (auto n : newPairToAdd) {
				if (std::find(nodeCheck.begin(), nodeCheck.end(), n.front())!=nodeCheck.end() &&
					 std::find(nodeCheck.begin(), nodeCheck.end(), n.back())!=nodeCheck.end()) {
					/*std::cout << "DEJA DANS LA LISTE DES NOEUDS " << n.front() << " " << n.back() << std::endl;
					std::cout<<"Pair check :" <<std::endl;
					for (auto i : nodeCheck){
						std::cout<<i<<std::endl;
					}*/
				}
				else {
					/*std::cout << "AJOUT : " << n.front() << " " << n.back() << std::endl;
					std::cout<<"Pair check :" <<std::endl;
					for (auto i : nodeCheck){
						std::cout<<i<<std::endl;
					}*/
					//std::cout<<"Element ADD : "<<n.front()<<" "<<n.back()<<std::endl;
					listPairToCheck.push_back(n);
				}
			}
			//Old method to remove the element
			//listPairToCheck.erase(std::remove(listPairToCheck.begin(), listPairToCheck.end(), p), listPairToCheck.end());

			listPairToCheck.erase(std::find(listPairToCheck.begin(),listPairToCheck.end(),p));

		}
	}
	/*
	for (auto n : listAllPairNodesChain){
		std::cout<<"\nUne Pair : "<<std::endl;
		for (auto i : n){
			std::cout<<i<<std::endl;
		}
	}*/
	return listAllPairNodesChain;
}

std::vector<Node> Tools::nodesAroundANode(const int i1)
{
	std::vector<Node> listNodesAround={};
	std::vector<Face> facesAround = Tools::getListFacesOfNode(i1);

	for (auto f : facesAround){
		auto listNodesFace = Tools::getListNodesOfFace(f.id());
		for (auto n : listNodesFace){
			if(Tools::checkExistEdge(i1,n.id(),f.id())){
				if(std::find(listNodesAround.begin(),listNodesAround.end(),n)!=listNodesAround.end()){
				}
				else{
					listNodesAround.push_back(n);
				}
			}
		}
	}
	/*std::cout<<"List Node Around"<<std::endl;
	for (auto n : listNodesAround){
		std::cout<<n<<std::endl;
	}*/

	return listNodesAround;
}

std::vector<Node> Tools::getOppositeNodes(const int i1,const int i2)
{
	std::vector<Node> oppositeNodes={};
	std::vector<Node> listAroundNode=Tools::nodesAroundANode(i1);
	auto opositeEdge = Tools::getOtherNodes(i1,i2);

	for (auto n : listAroundNode){
		for(auto v : opositeEdge){
			for(auto i : v){
				if (n.id()== i){
					if(std::find(oppositeNodes.begin(),oppositeNodes.end(),n)!=oppositeNodes.end()){
					}
					else{
						oppositeNodes.push_back(n);
					}
				}
			}
		}
	}
	/*std::cout<<"List Node Opposite"<<std::endl;
	for (auto n : oppositeNodes){
		std::cout<<n<<std::endl;
	}*/
	return oppositeNodes;
}

std::vector<int> Tools::getListFirstNodesChain(const int i1, const int i2)
{
	std::vector<int> listIdFirstNodes = {i1};
	auto listPairNodes = Tools::getAllNodesChain(i1,i2);
	auto nodesToCheck = Tools::getOppositeNodes(i1,i2);

	while(!nodesToCheck.empty()){
		auto listNodeParcours = nodesToCheck;
		for(auto n : listNodeParcours){
			listIdFirstNodes.push_back(n.id());
			for(auto v : listPairNodes){
				if(v.front() == n.id()){
					auto newNodeToAdd=Tools::getOppositeNodes(n.id(),v.back());
					for(auto add : newNodeToAdd){
						if(std::find(listIdFirstNodes.begin(),listIdFirstNodes.end(),add.id())!=listIdFirstNodes.end()){
						}
						else{
							nodesToCheck.push_back(add);
						}
					}
				}
				if(v.back()==n.id()){
					auto newNodeToAdd=Tools::getOppositeNodes(n.id(),v.front());
					for(auto add : newNodeToAdd){
						if(std::find(listIdFirstNodes.begin(),listIdFirstNodes.end(),add.id())!=listIdFirstNodes.end()){
						}
						else{
							nodesToCheck.push_back(add);
						}
					}
				}
			}
			//supp variable check
			nodesToCheck.erase(std::find(nodesToCheck.begin(),nodesToCheck.end(),n));
		}
	}

	/*std::cout<<"List Nodes First"<<std::endl;
	for (auto n : listIdFirstNodes){
		std::cout<<n<<std::endl;
	}*/

	return listIdFirstNodes;
}

std::vector<int> Tools::getListSecondNodesChain(const int i1, const int i2)
{
	auto listSecondNodesChain=Tools::getListFirstNodesChain(i2,i1);
	/*std::cout<<"List Nodes First"<<std::endl;
	for (auto n : listSecondNodesChain){
		std::cout<<n<<std::endl;
	}*/
	return listSecondNodesChain;
}

Node Tools::createMiddleNode(Node i1,Node i2)
{
	double newX = (i1.X()+i2.X())/2;
	double newY=(i1.Y()+i2.Y())/2;
	double newZ=(i1.Z()+i2.Z())/2;

	Node newNode = g_grid.m_mesh->newNode(math::Point(newX,newY,newZ));

	//std::cout<<"Le new node : " <<newNode<<std::endl;

	return newNode;
}

void Tools::createAllMiddlePoint(Node i1, Node i2)
{
	auto listPair = Tools::getAllNodesChain(i1.id(),i2.id());
	/*std::cout<<"List pair Nodes :"<<std::endl;
	for (auto n : listPair){
		for (auto i : n){
			std::cout<<i<<std::endl;
		}
	}*/
	for(auto p : listPair){
		Node nodeToAdd = Tools::createMiddleNode(g_grid.m_mesh->get<Node>(p.front()),g_grid.m_mesh->get<Node>(p.back()));
		g_grid.m_mesh->newNode(math::Point(nodeToAdd.X(),nodeToAdd.Y(),nodeToAdd.Z()));
		//std::cout<<"Le nouveau NOEUD :"<<nodeToAdd<<std::endl;
	}
}

std::vector<std::vector<int>> Tools::getPairNodesFace(const int i1, const int i2, const int faceID)
{
	std::vector<std::vector<int>> pairNodesFace = {};
	auto allPairNodes = Tools::getAllNodesChain(i1,i2);
	auto listNodes = Tools::getListNodesOfFace(faceID);

	/*for(auto p: allPairNodes){
		std::cout<<"Une pair mon gars :"<<std::endl;
		for(auto i : p){
			std::cout<<i<<std::endl;
		}
	}*/

	for (auto p : allPairNodes){
		if (std::find(listNodes.begin(),listNodes.end(),g_grid.m_mesh->get<Node>(p.front()))!=listNodes.end() ||
			 std::find(listNodes.begin(),listNodes.end(),g_grid.m_mesh->get<Node>(p.back()))!=listNodes.end()){
			pairNodesFace.push_back(p);
		}
	}

	/*
	std::cout<<"List pair Nodes of Face :"<<std::endl;
	for (auto n : pairNodesFace){
		for (auto i : n){
			std::cout<<i<<std::endl;
		}
	}*/
	if(pairNodesFace.empty()){
		std::cout<<"VIDE"<<std::endl;
	}

	return pairNodesFace;
}
std::vector<Node> Tools::getBoundaryEdge(Node firstNodeId, Node secondNodeId)
{
	auto listAllPair = Tools::getAllNodesChain(firstNodeId.id(),secondNodeId.id());
	std::vector<std::vector<Node>> boundaryEdges = {};
	std::vector<Node> listNodes = {};
	std::vector<Node> boundaryEdge = {};
	Node node;
	for (auto i : listAllPair) {
		auto faceCommon=Tools::getIdOneCommonFace(i.front(),i.back());
		auto facesCommon = Tools::getFacesCommon(i.front(),i.back());
		if(facesCommon.front().id()==faceCommon &&facesCommon.back().id()==faceCommon){
			std::vector<Node> edge ={};
			edge.push_back(g_grid.m_mesh->get<Node>(i.front()));
			edge.push_back(g_grid.m_mesh->get<Node>(i.back()));
			listNodes.push_back(g_grid.m_mesh->get<Node>(i.front()));
			listNodes.push_back(g_grid.m_mesh->get<Node>(i.back()));
			boundaryEdges.push_back(edge);
			/*std::cout<<"Edge Boundary :"<<std::endl;
			for (auto b : boundaryEdge) {
				std::cout<<b<<std::endl;
			}*/

		}
	}
	for (int i = 0; i <listNodes.size()-1;i++){
		for (int j = 0; j <listNodes.size()-1;j++){
			if(listNodes[j].X()>listNodes[j+1].X()){
				std::swap(listNodes[j],listNodes[j+1]);
			}
		}
	}
	/*
	std::cout<<"TRIEE :"<<std::endl;
	for (auto j : listNodes){
		std::cout<<j<<std::endl;
		std::cout<<j.X()<<std::endl;
	}*/
	if(listNodes[0].Y()<=listNodes[1].Y()){
		node = listNodes[0];
	}
	else{
		node = listNodes[1];
	}

	//std::cout<<node<<std::endl;

	for (auto e: boundaryEdges){
		if(e.front().X()==node.X() && e.front().Y()==node.Y()){
			boundaryEdge.push_back(e.front());
			boundaryEdge.push_back(e.back());
		}
		else if (e.back().X()==node.X()&&e.back().Y()==node.Y()){
			boundaryEdge.push_back(e.back());
			boundaryEdge.push_back(e.front());
		}
	}

	/*boundaryEdge.push_back(boundaryEdges.front().front());
	boundaryEdge.push_back(boundaryEdges.front().back());
	std::cout<<"L'arete finale :"<<std::endl;
	for (auto i:boundaryEdge){
		std::cout<<i<<std::endl;
	}*/

	return boundaryEdge;
}



void Tools::createEdge(Node i1, Node i2)
{
	auto edgeBoundary = Tools::getBoundaryEdge(i1,i2);
	auto listFaces = Tools::getAllFacesChain(edgeBoundary.front().id(),edgeBoundary.back().id());
	auto listEdges = Tools::getAllNodesChain(edgeBoundary.front().id(),edgeBoundary.back().id());
	auto firstNodesList = Tools::getListFirstNodesChain(edgeBoundary.front().id(),edgeBoundary.back().id());
	auto secondNodesList = Tools::getListSecondNodesChain(edgeBoundary.front().id(),edgeBoundary.back().id());
	std::vector<int> previousEdge;

	bool vertical = Tools::checkVertical(edgeBoundary.front(),edgeBoundary.back());


	Node previousMiddleNode;

	for (auto e : listEdges){

		if(e.front()!=edgeBoundary.front().id()&& e.back() != edgeBoundary.back().id()){

			auto faceEdges = Tools::getFaceEdges(previousEdge.front(),previousEdge.back(),e.front(),e.back());

			auto pairNodesFace = Tools::getPairNodesFace(edgeBoundary.front().id(), edgeBoundary.back().id(), faceEdges.id());
			std::vector<int> firstPair = pairNodesFace.front();
			std::vector<int> secondPair = pairNodesFace.back();

			int valueActivateFace;


			Node nodeFirstPair1;
			Node nodeSecondPair1;

			Node nodeFirstPair2;
			Node nodeSecondPair2;

			for (auto l : firstNodesList) {
				for (auto i : firstPair) {
					if (i == l) {
						nodeSecondPair1 = g_grid.m_mesh->get<Node>(i);
					}
				}
				for (auto i : secondPair) {
					if (i == l) {
						nodeSecondPair2 = g_grid.m_mesh->get<Node>(i);
					}
				}
			}
			for (auto l : secondNodesList) {
				for (auto i : firstPair) {
					if (i == l) {
						nodeFirstPair1 = g_grid.m_mesh->get<Node>(i);
					}
				}
				for (auto i : secondPair) {
					if (i == l) {
						nodeFirstPair2 = g_grid.m_mesh->get<Node>(i);
					}
				}
			}
			if(vertical){
				//std::cout<<"Vertical If"<<std::endl;
				Node newMiddleNode = Tools::createMiddleNode(g_grid.m_mesh->get<Node>(secondPair.front()), g_grid.m_mesh->get<Node>(secondPair.back()));
				Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair2,nodeFirstPair1,  previousMiddleNode,newMiddleNode);
				g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(newMiddleNode.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(previousMiddleNode.id()).add(newFaceLeft);

				/*
				std::cout <<"Liste Noeuds :"<<std::endl;
				std::cout<<"Node First Pair 2 "<<nodeFirstPair2<<std::endl;
				std::cout<<"Node First Pair 1 "<<nodeFirstPair1<<std::endl;
				std::cout<<"Previous Middle Node "<<previousMiddleNode<<std::endl;
				std::cout<<"newMiddleNode "<<newMiddleNode<<std::endl;
				std::cout<<"Face Left : "<< newFaceLeft<<std::endl;
				for (auto n : Tools::getListNodesOfFace(newFaceLeft.id())){
					std::cout<<n<<std::endl;
				}*/


				Face newFaceRight = g_grid.m_mesh->newQuad( newMiddleNode,previousMiddleNode,nodeSecondPair1,  nodeSecondPair2);
				g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(newMiddleNode.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(previousMiddleNode.id()).add(newFaceRight);

				/*
				std::cout <<"Liste Noeuds :"<<std::endl;

				std::cout<<"newMiddleNode "<<newMiddleNode<<std::endl;
				std::cout<<"Previous Middle Node "<<previousMiddleNode<<std::endl;
				std::cout<<"Node Second Pair 1 "<<nodeSecondPair1<<std::endl;
				std::cout<<"Node Second Pair 2 "<<nodeSecondPair2<<std::endl;
				std::cout<<"Face Right : "<< newFaceRight<<std::endl;
				for (auto n : Tools::getListNodesOfFace(newFaceRight.id())){
					std::cout<<n<<std::endl;
				}*/

				valueActivateFace=g_grid.m_mesh->getVariable<int,GMDS_FACE>("activate")->value(faceEdges.id());

				g_grid.activate->set(newFaceLeft.id(),g_grid.getActivate(faceEdges)); //attribution val aux faces
				g_grid.activate->set(newFaceRight.id(),g_grid.getActivate(faceEdges)); //attribution val aux faces

				previousMiddleNode = newMiddleNode;
				previousEdge = {e.front(), e.back()};
			}
			else{
				//std::cout<<"Horizontal If"<<std::endl;
				Node newMiddleNode = Tools::createMiddleNode(g_grid.m_mesh->get<Node>(secondPair.front()), g_grid.m_mesh->get<Node>(secondPair.back()));
				Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair2,  newMiddleNode,previousMiddleNode, nodeFirstPair1);
				g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(newMiddleNode.id()).add(newFaceLeft);
				g_grid.m_mesh->get<Node>(previousMiddleNode.id()).add(newFaceLeft);
				/*
				std::cout<<"Face Left : "<< newFaceLeft<<std::endl;
				for (auto n : Tools::getListNodesOfFace(newFaceLeft.id())){
					std::cout<<n<<std::endl;
				}*/


				Face newFaceRight = g_grid.m_mesh->newQuad( newMiddleNode,  nodeSecondPair2,nodeSecondPair1,previousMiddleNode);
				g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(newMiddleNode.id()).add(newFaceRight);
				g_grid.m_mesh->get<Node>(previousMiddleNode.id()).add(newFaceRight);
				/*
				std::cout<<"Face Right : "<< newFaceRight<<std::endl;
				for (auto n : Tools::getListNodesOfFace(newFaceRight.id())){
					std::cout<<n<<std::endl;
				}*/

				valueActivateFace=g_grid.m_mesh->getVariable<int,GMDS_FACE>("activate")->value(faceEdges.id());

				g_grid.activate->set(newFaceLeft.id(),g_grid.getActivate(faceEdges)); //attribution val aux faces
				g_grid.activate->set(newFaceRight.id(),g_grid.getActivate(faceEdges)); //attribution val aux faces

				previousMiddleNode = newMiddleNode;
				previousEdge = {e.front(), e.back()};

			}
		}
		else{
			previousEdge={e.front(),e.back()};
			previousMiddleNode=Tools::createMiddleNode(g_grid.m_mesh->get<Node>(e.front()),g_grid.m_mesh->get<Node>(e.back()));
		}
	}
	for (auto lFace : listFaces){
		auto listNodesOfFace = Tools::getListNodesOfFace(lFace.id());
		for (auto n : listNodesOfFace){
			g_grid.m_mesh->get<Node>(n.id()).remove(lFace);
		}
		g_grid.m_mesh->deleteFace(lFace.id());
	}

}
Face Tools::getFaceEdges(const int i1, const int i2, const int i3, const int i4)
{
	Face face;
	auto face1Edge = Tools::getFacesCommon(i1,i2);
	auto face2Edge = Tools::getFacesCommon(i3,i4);

	for (auto f1 : face1Edge){
		for (auto f2 : face2Edge){
			if (f1 == f2){
				return f1;
			}
		}
	}
}

bool
Tools::checkVertical(Node i1, Node i2)
{
	auto otherEdge = Tools::getOtherNodes(i1.id(),i2.id());
	auto faceEdge = g_grid.m_mesh->get<Face>(Tools::getIdOneCommonFace(i1.id(),i2.id()));
	auto boundaryEdge = Tools::getBoundaryEdge(i1,i2);
	auto listNodesCheck = Tools::getAllNodesSameWay(boundaryEdge.front(),boundaryEdge.back());

	for(auto i : listNodesCheck){
		//std::cout<<"LES NOEUDS VU : "<<i<<std::endl;
		if (i==1){
			//std::cout<<"VERTICAL"<<std::endl;
			return true ;
		}
	}
	//std::cout<<"HORIZONTAL"<<std::endl;
	return false;
}

std::vector<int> Tools::getAllNodesSameWay(Node i1, Node i2)
{

	auto otherNodes = Tools::getOtherNodes(i1.id(),i2.id());
	auto faceEdge = g_grid.m_mesh->get<Face>(Tools::getIdOneCommonFace(i1.id(),i2.id()));
	std::vector<int> edge={};
	if (Tools::checkExistEdge(i1.id(),otherNodes.front().front(),faceEdge.id())){
		edge.push_back(i1.id());
		edge.push_back(otherNodes.front().front());
	}
	else{
		edge.push_back(i1.id());
		edge.push_back(otherNodes.front().back());
	}
	std::vector<int> listNodes=Tools::getListFirstNodesChain(edge.front(),edge.back());
	/*std::cout<<"List Noeuds meme sens :"<<std::endl;
	for (auto i : listNodes){
		std::cout<<i<<std::endl;
	}*/
	return listNodes;
}

std::map<TCellID ,TCellID> Tools::getBoundaryNodes(Mesh *AMesh)
{

	gmds::MeshDoctor doc_imprint(AMesh);
	doc_imprint.updateUpwardConnectivity();
	doc_imprint.buildBoundaryCells();
	std::vector<Node> listBoundaryNodes={};

	Mesh boundaryMesh(MeshModel (DIM2|E|N|N2E|E2N));

	BoundaryExtractor2D boundary_extractor(AMesh,&boundaryMesh);
	std::map<TCellID ,TCellID > aNodeMap;
	std::map<TCellID ,TCellID > aEdgeMap;
	std::map<TCellID ,TCellID > aNodeMapInv;
	std::map<TCellID ,TCellID > aEdgeMapInv;
	boundary_extractor.setMappings(&aNodeMap,&aEdgeMap,&aNodeMapInv,&aEdgeMapInv);
	boundary_extractor.execute();

	return aNodeMap;
}
std::vector<Node> Tools::getVerticalEdge(Face AFace)
{
	std::vector<Node> verticalEdge;
	auto listNodesFace = Tools::getListNodesOfFace(AFace.id());
	/*std::cout<<"Liste Noeuds "<<std::endl;
	for (auto n : listNodesFace){
		std::cout<<"le noeud : "<<n<<std::endl;
	}*/
	for (auto i : listNodesFace){
		for (auto j : listNodesFace){
			if (i != j && Tools::checkExistEdge(i.id(),j.id(),AFace.id()) && Tools::checkVertical(i,j)){
				verticalEdge.push_back(i);
				verticalEdge.push_back(j);
				//std::cout<<"l'arete verticale : " << verticalEdge.front()<<" et " <<verticalEdge.back()<<std::endl;
				return verticalEdge;
			}
		}
	}
}
std::vector<Node> Tools::getHorizontalEdge(Face AFace)
{
	std::vector<Node> horizontalEdge;
	auto listNodesFace = Tools::getListNodesOfFace(AFace.id());
	for (auto i : listNodesFace){
		for (auto j : listNodesFace){
			if (i != j && Tools::checkExistEdge(i.id(),j.id(),AFace.id()) && !Tools::checkVertical(i,j)){
				horizontalEdge.push_back(i);
				horizontalEdge.push_back(j);
				//std::cout<<"l'arete horizontale : " << horizontalEdge.front()<<" et " <<horizontalEdge.back()<<std::endl;
				return horizontalEdge;
			}
		}
	}
	//std::cout<<"Arete select : "<<horizontalEdge.front()<<" \n"<<horizontalEdge.back()<<std::endl;
}

Node Tools::selectNodeMinRange(Face AFace)
{
	Node minNode;
	double minRange;
	bool minRangeNull = true;
	auto listNodes = Tools::getListNodesOfFace(AFace.id());
	auto  boundaryNodes = Tools::getBoundaryNodes(g_grid.meshTarget);
	for (auto b : boundaryNodes){
		for (auto n : listNodes) {
			double range = Tools::calcRangePoints(n,g_grid.meshTarget->get<Node>(b.first));
			if (minRangeNull == true ) {
				minNode = n;
				minRange = range;
				minRangeNull=false;
			}
			else if (range < minRange){
				minRange = range;
				minNode = n;

			}
		}

	}
	return minNode;

}

Node Tools::selectNodeMaxRange(Face AFace)
{
	Node maxNode;
	double maxRange;
	bool maxRangeNull = true;
	auto listNodes = Tools::getListNodesOfFace(AFace.id());
	auto  boundaryNodes = Tools::getBoundaryNodes(g_grid.meshTarget);
	for (auto b : boundaryNodes){
		for (auto n : listNodes) {
			double range = Tools::calcRangePoints(n,g_grid.meshTarget->get<Node>(b.first));
			if (maxRangeNull == true ) {
				maxNode = n;
				maxRange = range;
				maxRangeNull=false;
			}
			else if (range > maxRange){
				maxRange = range;
				maxNode = n;

			}
		}

	}
	return maxNode;

}

double Tools::calcRangePoints(Node node1, Node node2)
{
	double range = 0;
	range = sqrt(pow(node1.X()-node2.X(),2)+pow(node1.Z()-node2.Z(),2)+pow(node1.Y()-node2.Y(),2));
	return range;
}

void Tools::cloneMesh(const GridBuilderAround &originalGrid, GridBuilderAround &newGrid)
{
	if (originalGrid.getDim() == 3 or newGrid.getDim() == 3)
	{
		throw GMDSException("Dimension must be 2");
	}

	for (int nodeID:originalGrid.m_mesh->nodes())
	{
		newGrid.m_mesh->newNode(originalGrid.m_mesh->get<Node>(nodeID).point());
	}
	std::set<TCellID> invalidIndexes;
	for (int faceID = 0; faceID <= originalGrid.m_mesh->getMaxLocalID(2); faceID++)
	{
		if (originalGrid.m_mesh->has<Face>(faceID))
		{
			newGrid.m_mesh->newFace(originalGrid.m_mesh->get<Face>(faceID).get<Node>());
			newGrid.activate->set(faceID,originalGrid.m_mesh->getVariable<int,GMDS_FACE>("activate")->value(faceID));
			newGrid.volFrac->set(faceID,originalGrid.m_mesh->getVariable<double,GMDS_FACE>("volFrac")->value(faceID));
			auto listNodeFaceOriginalGrid=Tools::getListNodesOfFace(faceID);
			for (int i =0;i<4;i++){
				auto actualNode = listNodeFaceOriginalGrid[i];
				auto actualFace = newGrid.m_mesh->get<Face>(faceID);
				newGrid.m_mesh->get<Node>(actualNode.id()).add(actualFace);
			}
		}
		else
		{
			newGrid.m_mesh->newFace({0, 0, 0, 0});
			invalidIndexes.insert(faceID);
		}
	}
	for (int faceID:invalidIndexes)
	{
		newGrid.m_mesh->deleteFace(faceID);
	}
}
//========================================================================================
//OLD FUNCTION

void Tools::joinFaceToNodes(Node i1, Node i2)
{
	auto boundaryEdge = Tools::getBoundaryEdge(i1, i2);
	auto listFaces = Tools::getAllFacesChain(boundaryEdge.front().id(), boundaryEdge.back().id());
	auto listPairNodes = Tools::getAllNodesChain(boundaryEdge.front().id(), boundaryEdge.back().id());
	auto listFirstNodes = Tools::getListFirstNodesChain(boundaryEdge.front().id(), boundaryEdge.back().id());
	auto listSecondNodes = Tools::getListSecondNodesChain(boundaryEdge.front().id(), boundaryEdge.back().id());

	std::vector<Node> nodesCreate;
	std::vector<Node> nodesDouble;


	//Node oldNewNodePair1 = g_grid.m_mesh->newNode(10000000,10000000,10000000);
	//Node oldNewNodePair2 = g_grid.m_mesh->newNode(10000000,10000000,10000000);

	for (auto f : listFaces) {
		auto pairNodesFace = Tools::getPairNodesFace(i1.id(), i2.id(), f.id());
		std::vector<int> firstPair = pairNodesFace.front();
		std::vector<int> secondPair = pairNodesFace.back();

		auto newNodePair1 = Tools::createMiddleNode(g_grid.m_mesh->get<Node>(firstPair.front()), g_grid.m_mesh->get<Node>(firstPair.back()));
		auto newNodePair2 = Tools::createMiddleNode(g_grid.m_mesh->get<Node>(secondPair.front()), g_grid.m_mesh->get<Node>(secondPair.back()));


		Node nodeFirstPair1;
		Node nodeSecondPair1;

		Node nodeFirstPair2;
		Node nodeSecondPair2;

		for (auto l : listFirstNodes) {
			for (auto i : firstPair) {
				if (i == l) {
					nodeSecondPair1 = g_grid.m_mesh->get<Node>(i);
				}
			}
			for (auto i : secondPair) {
				if (i == l) {
					nodeSecondPair2 = g_grid.m_mesh->get<Node>(i);
				}
			}
		}
		for (auto l : listSecondNodes) {
			for (auto i : firstPair) {
				if (i == l) {
					nodeFirstPair1 = g_grid.m_mesh->get<Node>(i);
				}
			}
			for (auto i : secondPair) {
				if (i == l) {
					nodeFirstPair2 = g_grid.m_mesh->get<Node>(i);
				}
			}
		}

		/*if (newNodePair1.X() == oldNewNodePair1.X() && newNodePair1.Y() == oldNewNodePair1.Y() && newNodePair1.Z() == oldNewNodePair1.Z()) {
		   std::cout << "DANS LE IF 1 :" << std::endl;
		   Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, oldNewNodePair1, newNodePair2, nodeFirstPair2);
		   g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		   g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		   g_grid.m_mesh->get<Node>(oldNewNodePair1.id()).add(newFaceLeft);
		   g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceLeft);

		   Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, oldNewNodePair1, newNodePair2, nodeSecondPair2);
		   g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		   g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		   g_grid.m_mesh->get<Node>(oldNewNodePair1.id()).add(newFaceRight);
		   g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceRight);

		   oldNewNodePair1=newNodePair1;
		   oldNewNodePair2=newNodePair2;
		}
		else {
		   if (newNodePair1.X() == oldNewNodePair2.X() && newNodePair1.Y() == oldNewNodePair2.Y() && newNodePair1.Z() == oldNewNodePair2.Z()) {
		      std::cout << "DANS LE IF 2 :" << std::endl;
		      Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, oldNewNodePair2, newNodePair2, nodeFirstPair2);
		      g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		      g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		      g_grid.m_mesh->get<Node>(oldNewNodePair2.id()).add(newFaceLeft);
		      g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceLeft);

		      Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, oldNewNodePair2, newNodePair2, nodeSecondPair2);
		      g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		      g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		      g_grid.m_mesh->get<Node>(oldNewNodePair2.id()).add(newFaceRight);
		      g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceRight);
		      oldNewNodePair1=newNodePair1;
		      oldNewNodePair2=newNodePair2;
		   }
		   else {
		      if (newNodePair2.X() == oldNewNodePair1.X() && newNodePair2.Y() == oldNewNodePair1.Y() && newNodePair2.Z() == oldNewNodePair1.Z()) {
		         std::cout << "DANS LE IF 3 :" << std::endl;
		         Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, newNodePair1, oldNewNodePair1, nodeFirstPair2);
		         g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		         g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		         g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceLeft);
		         g_grid.m_mesh->get<Node>(oldNewNodePair1.id()).add(newFaceLeft);

		         Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, newNodePair1, oldNewNodePair1, nodeSecondPair2);
		         g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		         g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		         g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceRight);
		         g_grid.m_mesh->get<Node>(oldNewNodePair1.id()).add(newFaceRight);
		         oldNewNodePair1=newNodePair1;
		         oldNewNodePair2=newNodePair2;
		      }
		      else {
		         if (newNodePair2.X() == oldNewNodePair2.X() && newNodePair2.Y() == oldNewNodePair2.Y() && newNodePair2.Z() == oldNewNodePair2.Z()) {
		            std::cout << "DANS LE IF 4 :" << std::endl;
		            Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, newNodePair1, oldNewNodePair2, nodeFirstPair2);
		            g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(oldNewNodePair2.id()).add(newFaceLeft);

		            Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, newNodePair1, oldNewNodePair2, nodeSecondPair2);
		            g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(oldNewNodePair2.id()).add(newFaceRight);
		            oldNewNodePair1=newNodePair1;
		            oldNewNodePair2=newNodePair2;
		         }
		         else {
		            std::cout << "DANS LE ELSE FINAL :" << std::endl;
		            Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, newNodePair1, newNodePair2, nodeFirstPair2);
		            g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceLeft);
		            g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceLeft);

		            Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, newNodePair1, newNodePair2, nodeSecondPair2);
		            g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceRight);
		            g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceRight);

		            oldNewNodePair1=newNodePair1;
		            oldNewNodePair2=newNodePair2;

		         }
		      }
		   }
		}
		*/
		/*if(!nodesCreate.empty()) {
		   for (auto w : nodesCreate) {
		      if (w.X() == newNodePair1.X() && w.Y() == newNodePair1.Y() && w.Z() == newNodePair1.Z()) {
		         nodesDouble.push_back(newNodePair1);
		      }
		      if (w.X() == newNodePair2.X() && w.Y() == newNodePair2.Y() && w.Z() == newNodePair2.Z()) {
		         nodesDouble.push_back(newNodePair2);
		      }
		      if (w.X() != newNodePair1.X() && w.Y() != newNodePair1.Y() && w.Z() != newNodePair1.Z()) {
		      }
		   }
		}*/
		std::cout << "DANS LE ELSE FINAL :" << std::endl;
		Face newFaceLeft = g_grid.m_mesh->newQuad(nodeFirstPair1, newNodePair1, newNodePair2, nodeFirstPair2);
		g_grid.m_mesh->get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
		g_grid.m_mesh->get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
		g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceLeft);
		g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceLeft);

		Face newFaceRight = g_grid.m_mesh->newQuad(nodeSecondPair1, newNodePair1, newNodePair2, nodeSecondPair2);
		g_grid.m_mesh->get<Node>(nodeSecondPair1.id()).add(newFaceRight);
		g_grid.m_mesh->get<Node>(nodeSecondPair2.id()).add(newFaceRight);
		g_grid.m_mesh->get<Node>(newNodePair1.id()).add(newFaceRight);
		g_grid.m_mesh->get<Node>(newNodePair2.id()).add(newFaceRight);

		for (auto n : nodesCreate) {
			if (n.X() == newNodePair1.X() && n.Y() == newNodePair1.Y() && n.Z() == newNodePair1.Z()) {
				nodesDouble.push_back(newNodePair1);
			}
			if (n.X() == newNodePair2.X() && n.Y() == newNodePair2.Y() && n.Z() == newNodePair2.Z()) {
				nodesDouble.push_back(newNodePair2);
			}
		}

		nodesCreate.push_back(newNodePair1);
		nodesCreate.push_back(newNodePair2);

	}
	std::cout<<"NODES CREE EN DOUBLE :"<<std::endl;
	for (auto n : nodesDouble){
		std::cout<<n<<std::endl;
		//g_grid.m_mesh->deleteNode(n.id());
	}
	for (auto lFace : listFaces){
		g_grid.m_mesh->deleteFace(lFace.id());
	}
}


void Tools::intersectionTargetWithGrid(Mesh *AMesh, const Mesh *AImprintMesh, gmds::Variable<double> *AVolFrac)
{
	// check validity of the inputs
	bool valid_input = true;
	std::string msg("volfraccomputation_2d ");

	gmds::MeshModel model = AMesh->getModel();
	if (!model.has(gmds::F) || !model.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AMesh->getNbFaces() != AMesh->getNbQuadrilaterals()) {
		msg += std::string("AMesh should have only quads");
		valid_input = false;
	}
	gmds::MeshModel modelImprint = AImprintMesh->getModel();
	if (!modelImprint.has(gmds::F) || !modelImprint.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AImprintMesh->getNbFaces() != AImprintMesh->getNbTriangles()) {
		msg += std::string("AImprintMesh should have only triangles");
		valid_input = false;
	}

	// check mesh orientation
	// TODO
	for(auto f_id: AMesh->faces()) {
		gmds::Face f = AMesh->get<Face>(f_id);
		std::vector<gmds::Node> n = f.get<gmds::Node>();
		AVolFrac->set(f_id,0);

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		//std::cout<<"SJ Face "<<f_id<<"("<<n[0].id()<<", "<<n[1].id()<<", "<<n[2].id()<<", "<<n[3].id()<<": "<<sj<<std::endl;
		if(sj < 0) {
			msg += std::string("AMesh has a bad cell.");
			valid_input = false;
			break;
		}
	}
/*
	for(auto f_id_imprint: AImprintMesh->faces()) {
		gmds::Face f = AImprintMesh->get<Face>(f_id_imprint);
		std::vector<gmds::Node> n = f.get<gmds::Node>();
		gmds::math::Triangle tri(n[0].point(), n[1].point(), n[2].point());
		double sj = tri.computeScaledJacobian2D();
		if(sj < 0) {
			msg += std::string("AImprintMesh has a bad cell.");
			valid_input = false;
			break;
		}
	}*/

	for (auto f_id: AMesh->faces()){

		std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

		std::cout<<"Face Grid select"<<f_id<<std::endl;

		std::cout<<"\"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\" "<<std::endl;


		gmds::Face f = AMesh->get<Face>(f_id);

		r2d_rvec2 vertices_Quad[4];
		std::vector<gmds::Node> n_quad = f.get<gmds::Node>();
		gmds::TCoord xyz[4][3];
		xyz[0][0] = n_quad[0].X();
		xyz[0][1] = n_quad[0].Y();
		xyz[0][2] = n_quad[0].Z();
		xyz[1][0] = n_quad[1].X();
		xyz[1][1] = n_quad[1].Y();
		xyz[1][2] = n_quad[1].Z();
		xyz[2][0] = n_quad[2].X();
		xyz[2][1] = n_quad[2].Y();
		xyz[2][2] = n_quad[2].Z();
		xyz[3][0] = n_quad[3].X();
		xyz[3][1] = n_quad[3].Y();
		xyz[3][2] = n_quad[3].Z();
		for(auto tri_id: AImprintMesh->faces()) {
			gmds::Face tri = AImprintMesh->get<Face>(tri_id);
			std::vector<gmds::Node> n_tri = tri.get<gmds::Node>();

			r2d_rvec2 vertices[3];
			vertices[0].x = n_tri[0].X();
			vertices[0].y = n_tri[0].Y();
			vertices[1].x = n_tri[2].X();
			vertices[1].y = n_tri[2].Y();
			vertices[2].x = n_tri[1].X();
			vertices[2].y = n_tri[1].Y();
			r2d_int numverts = 3;
			r2d_plane planes[3];
			r2d_poly_faces_from_verts(planes, vertices, numverts);

			gmds::math::Point tri_pt0(vertices[0].x,vertices[0].y);
			gmds::math::Point tri_pt1(vertices[1].x,vertices[1].y);
			gmds::math::Point tri_pt2(vertices[2].x,vertices[2].y);
			gmds::math::Triangle t(tri_pt0,tri_pt1,tri_pt2);




			r2d_int numverts_Quad = 4;
			r2d_plane planes_Quad[4];
			r2d_poly_faces_from_verts(planes_Quad, vertices_Quad, numverts_Quad);

			gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
			gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
			gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
			gmds::math::Point pt3(xyz[3][0], xyz[3][1], xyz[3][2]);
			gmds::math::Quadrilateral q(pt0, pt1, pt2, pt3);

			const double volsurf = q.area();
			const double IoU = t.area();


			r2d_poly poly;
			r2d_rvec2 verts[4];

			r2d_poly tetra;
			r2d_rvec2 verts_tri[3];


			r2d_init_poly(&tetra,vertices,3);

			r2d_clip(&tetra,planes_Quad,4);
			r2d_int POLY_ORDER = 2;
			r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
			r2d_reduce(&tetra, om, POLY_ORDER);

			double vf_tri = AVolFrac->value(tri_id);

			std::cout<<"Id Face : "<<tri_id<<std::endl;
			std::cout<<"valeur om : "<<om<<std::endl;
			std::cout<<"valeur om[0] : "<<om[0]<<std::endl;
			std::cout<<"valeur IoU air : "<<IoU<<std::endl;
			std::cout<<"old value : "<<vf_tri<<std::endl;



			AVolFrac->set(tri_id, vf_tri+om[0]*IoU);

			std::cout<<"valeur volfrac : "<<AVolFrac->value(tri_id)<<std::endl;



		}
	}
}
void Tools::anotherVolFrac(Mesh *AMesh, const Mesh *AImprintMesh, gmds::Variable<double> *AVolFrac)
{
	// check validity of the inputs
	bool valid_input = true;
	std::string msg("volfraccomputation_2d ");

	gmds::MeshModel model = AMesh->getModel();
	if (!model.has(gmds::F) || !model.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}
	if(AMesh->getNbFaces() != AMesh->getNbQuadrilaterals()) {
		msg += std::string("AMesh should have only quads");
		valid_input = false;
	}
	gmds::MeshModel modelImprint = AImprintMesh->getModel();
	if (!modelImprint.has(gmds::F) || !modelImprint.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}
	if(AImprintMesh->getNbFaces() != AImprintMesh->getNbTriangles()) {
		msg += std::string("AImprintMesh should have only triangles");
		valid_input = false;
	}

	// check mesh orientation
	// TODO
	for(auto f_id: AMesh->faces()) {
		AVolFrac->set(f_id, 0);
		gmds::Face f = AMesh->get<Face>(f_id);
		std::vector<gmds::Node> n = f.get<gmds::Node>();

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		/*if(sj < 0) {
			msg += std::string("AMesh has a bad cell.");
			valid_input = false;
			break;
		}*/
	}

	for(auto tri_id: AImprintMesh->faces()) {
		AVolFrac->set(tri_id, 0);
	}

	//	for(auto f_id: AImprintMesh->faces()) {
	//		gmds::Face f = AImprintMesh->get<Face>(f_id);
	//		std::vector<gmds::Node> n = f.get<gmds::Node>();
	//
	//		gmds::math::Triangle tri(n[0].point(), n[1].point(), n[2].point());
	//		double sj = tri.computeScaledJacobian2D();
	//		if(sj < 0) {
	//			msg += std::string("AImprintMesh has a bad cell.");
	//			valid_input = false;
	//			break;
	//		}
	//	}

	if(!valid_input) {
		throw gmds::GMDSException(msg);
	}

	for(auto tri_id: AImprintMesh->faces()) {
		gmds::Face tri = AImprintMesh->get<Face>(tri_id);
		std::vector<gmds::Node> n_tri = tri.get<gmds::Node>();

		r2d_rvec2 vertices[3];
		vertices[0].x = n_tri[0].X();
		vertices[0].y = n_tri[0].Y();
		vertices[1].x = n_tri[2].X();
		vertices[1].y = n_tri[2].Y();
		vertices[2].x = n_tri[1].X();
		vertices[2].y = n_tri[1].Y();
		r2d_int numverts = 3;
		r2d_plane planes[3];
		r2d_poly_faces_from_verts(planes,  vertices, numverts);

		for(auto f_id: AMesh->faces()) {
			gmds::Face f = AMesh->get<Face>(f_id);
			if(getValueActivateFace(f.id())==1) {

				std::vector<gmds::Node> n_quad = f.get<gmds::Node>();
				gmds::TCoord xyz[4][3];
				xyz[0][0] = n_quad[0].X();
				xyz[0][1] = n_quad[0].Y();
				xyz[0][2] = n_quad[0].Z();
				xyz[1][0] = n_quad[1].X();
				xyz[1][1] = n_quad[1].Y();
				xyz[1][2] = n_quad[1].Z();
				xyz[2][0] = n_quad[2].X();
				xyz[2][1] = n_quad[2].Y();
				xyz[2][2] = n_quad[2].Z();
				xyz[3][0] = n_quad[3].X();
				xyz[3][1] = n_quad[3].Y();
				xyz[3][2] = n_quad[3].Z();

				gmds::math::Point pt0(vertices[0].x, vertices[0].y);
				gmds::math::Point pt1(vertices[1].x, vertices[1].y);
				gmds::math::Point pt2(vertices[2].x, vertices[2].y);
				gmds::math::Triangle t(pt0, pt1, pt2);

				const double volsurf = t.area();

				r2d_poly poly;
				r2d_rvec2 verts[4];

				verts[0].x = xyz[0][0];
				verts[0].y = xyz[0][1];
				verts[1].x = xyz[1][0];
				verts[1].y = xyz[1][1];
				verts[2].x = xyz[2][0];
				verts[2].y = xyz[2][1];
				verts[3].x = xyz[3][0];
				verts[3].y = xyz[3][1];

				r2d_init_poly(&poly, verts, 4);

				r2d_clip(&poly, planes, 3);
				r2d_int POLY_ORDER = 2;
				r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
				r2d_reduce(&poly, om, POLY_ORDER);

				double vf = AVolFrac->value(tri_id);
				/*
				         std::cout<<"Id Face : "<<tri_id<<std::endl;
				         std::cout<<"valeur om : "<<om<<std::endl;
				         std::cout<<"valeur om[0] : "<<om[0]<<std::endl;
				         std::cout<<"valeur IoU : "<<volsurf<<std::endl;
				         std::cout<<"old value : "<<vf<<std::endl;*/

				//			std::cout<<"f_id "<<f_id<<" volsurf "<<volsurf<<" vf "<<vf<<" om "<<om[0]/volsurf<<std::endl;
				// add but do not forget to divide by cell area
				AVolFrac->set(tri_id, vf + om[0] / volsurf);
			}
		}
	}
}
void Tools::saveMesh(Mesh* AMesh,std::string filename)
{
	IGMeshIOService ioService(AMesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write(filename);
}

int Tools::bestCutDirection(gmds::Face AFace)
{
	auto listNodes = getListNodesOfFace(AFace.id());
	Node p0 = listNodes[0];
	Node p1 = listNodes[1];
	Node p2 = listNodes[2];
	Node p3 = listNodes[3];

	Node p5= createMiddleNode(p0,p1);
	Node p6= createMiddleNode(p1,p2);
	Node p7= createMiddleNode(p2,p3);
	Node p8= createMiddleNode(p3,p0);

	Mesh m0(MeshModel(DIM2|F|N|F2N|N2F));
	Mesh m1(MeshModel(DIM2|F|N|F2N|N2F));
	Mesh m2(MeshModel(DIM2|F|N|F2N|N2F));
	Mesh m3(MeshModel(DIM2|F|N|F2N|N2F));

	//mesh 0
	Node realp0m0=m0.newNode(p0.X(),p0.Y(),p0.Z());
	Node realp1m0=m0.newNode(p1.X(),p1.Y(),p1.Z());
	Node realp2m0=m0.newNode(p2.X(),p2.Y(),p2.Z());
	Node realp3m0=m0.newNode(p3.X(),p3.Y(),p3.Z());
	Node realp5m0=m0.newNode(p5.X(),p5.Y(),p5.Z());
	Node realp6m0=m0.newNode(p6.X(),p6.Y(),p6.Z());
	Node realp7m0=m0.newNode(p7.X(),p7.Y(),p7.Z());
	Node realp8m0=m0.newNode(p8.X(),p8.Y(),p8.Z());

	//mesh 1
	Node realp0m1=m1.newNode(p0.X(),p0.Y(),p0.Z());
	Node realp1m1=m1.newNode(p1.X(),p1.Y(),p1.Z());
	Node realp2m1=m1.newNode(p2.X(),p2.Y(),p2.Z());
	Node realp3m1=m1.newNode(p3.X(),p3.Y(),p3.Z());
	Node realp5m1=m1.newNode(p5.X(),p5.Y(),p5.Z());
	Node realp6m1=m1.newNode(p6.X(),p6.Y(),p6.Z());
	Node realp7m1=m1.newNode(p7.X(),p7.Y(),p7.Z());
	Node realp8m1=m1.newNode(p8.X(),p8.Y(),p8.Z());

	//mesh 2
	Node realp0m2=m2.newNode(p0.X(),p0.Y(),p0.Z());
	Node realp1m2=m2.newNode(p1.X(),p1.Y(),p1.Z());
	Node realp2m2=m2.newNode(p2.X(),p2.Y(),p2.Z());
	Node realp3m2=m2.newNode(p3.X(),p3.Y(),p3.Z());
	Node realp5m2=m2.newNode(p5.X(),p5.Y(),p5.Z());
	Node realp6m2=m2.newNode(p6.X(),p6.Y(),p6.Z());
	Node realp7m2=m2.newNode(p7.X(),p7.Y(),p7.Z());
	Node realp8m2=m2.newNode(p8.X(),p8.Y(),p8.Z());

	//mesh 3
	Node realp0m3=m3.newNode(p0.X(),p0.Y(),p0.Z());
	Node realp1m3=m3.newNode(p1.X(),p1.Y(),p1.Z());
	Node realp2m3=m3.newNode(p2.X(),p2.Y(),p2.Z());
	Node realp3m3=m3.newNode(p3.X(),p3.Y(),p3.Z());
	Node realp5m3=m3.newNode(p5.X(),p5.Y(),p5.Z());
	Node realp6m3=m3.newNode(p6.X(),p6.Y(),p6.Z());
	Node realp7m3=m3.newNode(p7.X(),p7.Y(),p7.Z());
	Node realp8m3=m3.newNode(p8.X(),p8.Y(),p8.Z());



	Face fm0 = m0.newQuad(realp0m0,realp5m0,realp7m0,realp3m0);
	Face f1m0 = m0.newQuad(realp5m0,realp1m0,realp2m0,realp7m0);
	Face fm1 = m1.newQuad(realp5m1,realp1m1,realp2m1,realp7m1);
	Face f1m1 = m1.newQuad(realp0m1,realp5m1,realp7m1,realp3m1);
	Face fm2 = m2.newQuad(realp8m2,realp6m2,realp2m2,realp3m2);
	Face f1m2 = m2.newQuad(realp0m2,realp1m2,realp6m2,realp8m2);
	Face fm3 = m3.newQuad(realp0m3,realp1m3,realp6m3,realp8m3);
	Face f1m3 = m3.newQuad(realp8m3,realp6m3,realp2m3,realp3m3);


	Variable<double>* volFracM0 = m0.newVariable<double,GMDS_FACE>("volFrac");
	Variable<double>* volFracM1 = m1.newVariable<double,GMDS_FACE>("volFrac");
	Variable<double>* volFracM2 = m2.newVariable<double,GMDS_FACE>("volFrac");
	Variable<double>* volFracM3 = m3.newVariable<double,GMDS_FACE>("volFrac");


	volfraccomputation_2d(&m0,g_grid.meshTarget,m0.getVariable<double,GMDS_FACE>("volFrac"));
	double IoUm0 = m0.getVariable<double,GMDS_FACE>("volFrac")->value(fm0.id());
	//std::cout<<"Valeur IoUm0 : "<<IoUm0<<std::endl;

	volfraccomputation_2d(&m1,g_grid.meshTarget,volFracM1);
	double IoUm1 = m1.getVariable<double,GMDS_FACE>("volFrac")->value(fm1.id());
	//std::cout<<"Valeur IoUm1 : "<<IoUm1<<std::endl;

	volfraccomputation_2d(&m2,g_grid.meshTarget,volFracM2);
	double IoUm2 = m2.getVariable<double,GMDS_FACE>("volFrac")->value(fm2.id());
	//std::cout<<"Valeur IoUm2 : "<<IoUm2<<std::endl;

	volfraccomputation_2d(&m3,g_grid.meshTarget,volFracM3);
	double IoUm3 = m3.getVariable<double,GMDS_FACE>("volFrac")->value(fm3.id());
	//std::cout<<"Valeur IoUm3 : "<<IoUm3<<std::endl;

	std::vector<double> listIoU{IoUm0,IoUm1,IoUm2,IoUm3};
	double maxIoU = 0;
	int cutSelect;


	for (int i = 0;i<listIoU.size();i++){
		//std::cout<<"IoU de : "<<i<<" est egal a : "<<listIoU[i]<<std::endl;
		if (listIoU[i]>maxIoU){
			maxIoU=listIoU[i];
			if(i==0 || i==1){
				cutSelect=1;
			}
			else{
				cutSelect=0;
			}
		}
	}
	//std::cout<<"Cut select : "<<cutSelect<<std::endl;
/*
	saveMesh(&m0,"m0save.vtk");
	saveMesh(&m1,"m1save.vtk");
	saveMesh(&m2,"m2save.vtk");
	saveMesh(&m3,"m3save.vtk");*/

	return cutSelect;

}

bool Tools::isValid(gmds::Face& AFace){
	std::vector<gmds::Node> nodes = AFace.get<Node>();
	Node nlast = nodes[nodes.size()-1];
	Node next = nodes[1];

	math::Point pc=nodes[0].point();
	math::Point pp=nlast.point();
	math::Point pn=next.point();
	math::Vector3d vprev(pc,pp);
	math::Vector3d vnext(pc,pn);
	math::Vector3d v0= vprev.cross(vnext);

	for(auto i=1;i<nodes.size();i++){

		Node prev= nodes[i-1];
		Node next = nodes[(i+1)%nodes.size()];

		pc=nodes[i].point();
		pp=prev.point();
		pn=next.point();
		math::Vector3d vprev(pc,pp);
		math::Vector3d vnext(pc,pn);
		math::Vector3d v= vprev.cross(vnext);
		if(v.dot(v0)<0)
			return false;
	}

	return true;
}

bool Tools::commonNodesFaces(gmds::Face AFace0, gmds::Face AFace1)
{
	auto listNodesFace0 = getListNodesOfFace(AFace0.id());
	auto listNodesFace1 = getListNodesOfFace(AFace1.id());
	int nbNodesCommon=0;

	for (auto nFace0 : listNodesFace0){
		for (auto nFace1 : listNodesFace1){
			if(nFace0.id()==nFace1.id()){
				nbNodesCommon++;
			}
		}
	}

	if(nbNodesCommon==2){
		std::cout<<"CEST VRAI"<<std::endl;
		return true;
	}
	else{
		std::cout<<"CEST FAUX"<<std::endl;
		return false;}
}
