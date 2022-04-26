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
	Face f = g_grid.m_mesh.get<Face>(faceID);
	std::vector<Node>list_nodes = f.get<Node>();
	/*for (auto n : list_nodes){
		std::cout<<"Face :"<< faceID <<" avec les noeuds :"<< n<<std::endl;
	}*/
	return list_nodes;
}

std::vector<Face> Tools::getListFacesOfNode(const int nodeID)
{
	Node n =g_grid.m_mesh.get<Node>(nodeID);
	std::vector<Face> list_faces = n.get<Face>();
	/*for (auto n : list_faces){
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
	std::vector<Face> listAllFaces = Tools::getFacesCommon(i1,i2);
	std::vector<std::vector<int>> listAllPairNodes = Tools::getOtherNodes(i1,i2);
	std::vector<int> nodeCheck = {i1,i2};
	auto listPairToCheck = listAllPairNodes;

	while(!listPairToCheck.empty()){
		auto listParcour = listPairToCheck;
		for (auto p : listParcour){
			/*std::cout<<"LES ID QUI SONT CHECK"<<std::endl;
			for (auto n : p){
				std::cout<<n<<std::endl;
			}*/
			if (*find(nodeCheck.begin(), nodeCheck.end(), p.front()) == p.front() ||
				 *find(nodeCheck.begin(), nodeCheck.end(), p.back()) == p.back()) {
				//std::cout<<"Noeud DEJA Check"<<std::endl;

			}
			else{
				nodeCheck.push_back(p.front());
				nodeCheck.push_back(p.back());

			}
			/*std::cout<<"Noeuds Check :"<<std::endl;
			for (auto j : nodeCheck){
				std::cout<<j<<std::endl;
			}*/
			std::vector<std::vector<int>> otherNodes = Tools::getOtherNodes(p.front(),p.back());


			/*for(auto f : Tools::getFacesCommon(p.front(),p.back())){
				listAllFaces.push_back(f);
				std::cout<<"Face ADD : "<<f<<std::endl;
			}*/
			for(auto n : otherNodes){
				if (*find(nodeCheck.begin(), nodeCheck.end(), n.front()) == n.front() ||
					 *find(nodeCheck.begin(), nodeCheck.end(), n.back()) == n.back()) {
					//std::cout<<"Other nodes deja check"<<std::endl;

				}
				else{
					listPairToCheck.push_back(n);
				}
			}
			auto newFaces = Tools::getFacesCommon(p.front(),p.back());
			for (auto f : newFaces){

				if (*find(listAllFaces.begin(), listAllFaces.end(),f) == f){
					//std::cout<<"Other nodes deja check"<<std::endl;
				}
				else{
					listAllFaces.push_back(f);
				}
			}

			//Old method to remove an element
			// listPairToCheck.erase(std::remove(listPairToCheck.begin(), listPairToCheck.end(), p), listPairToCheck.end());

			listPairToCheck.erase(std::find(listPairToCheck.begin(),listPairToCheck.end(),p));

			for (auto n : listPairToCheck){
				for (auto j : n){
					//std::cout<<"LIST EN COURS DE SUPP : "<<j<<std::endl;
				}
			}

		}
	}
	/*std::cout<<"LA LISTE DES FACES "<<std::endl;
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
		for (auto p : listParcour) {
			listAllPairNodesChain.push_back(p);
			pairNodes.push_back(p.front());
			pairNodes.push_back(p.back());
			auto newPairToAdd = Tools::getOtherNodes(p.front(), p.back());
			for (auto n : newPairToAdd) {
				if ((*find(pairNodes.begin(), pairNodes.end(), n.front()) == n.front()) ||
				    (*find(pairNodes.begin(), pairNodes.end(), n.back()) == n.back())) {
					std::cout << "DEJA DANS LA LISTE DES NOEUDS " << n.front() << " " << n.back() << std::endl;
					std::cout<<"Pair check :" <<std::endl;
					for (auto i : pairNodes){
						std::cout<<i<<std::endl;
					}
				}
				else {

					//std::cout<<"Element ADD : "<<n.front()<<" "<<n.back()<<std::endl;
					listPairToCheck.push_back(n);
				}
			}
			//Old method to remove the element
			//listPairToCheck.erase(std::remove(listPairToCheck.begin(), listPairToCheck.end(), p), listPairToCheck.end());

			listPairToCheck.erase(std::find(listPairToCheck.begin(),listPairToCheck.end(),p));

		}
	}

	for (auto n : listAllPairNodesChain){
		std::cout<<"\nUne Pair : "<<std::endl;
		for (auto i : n){
			std::cout<<i<<std::endl;
		}
	}
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

	Node newNode = g_grid.m_mesh.newNode(math::Point(newX,newY,newZ));

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
		Node nodeToAdd = Tools::createMiddleNode(g_grid.m_mesh.get<Node>(p.front()),g_grid.m_mesh.get<Node>(p.back()));
		g_grid.m_mesh.newNode(math::Point(nodeToAdd.X(),nodeToAdd.Y(),nodeToAdd.Z()));
		//std::cout<<"Le nouveau NOEUD :"<<nodeToAdd<<std::endl;
	}
}

std::vector<std::vector<int>> Tools::getPairNodesFace(const int i1, const int i2, const int faceID)
{
	std::vector<std::vector<int>> pairNodesFace = {};
	auto allPairNodes = Tools::getAllNodesChain(i1,i2);
	auto listNodes = Tools::getListNodesOfFace(faceID);

	for(auto p: allPairNodes){
		std::cout<<"Une pair :"<<std::endl;
		for(auto i : p){
			std::cout<<i<<std::endl;
		}
	}

	for (auto p : allPairNodes){
		if (std::find(listNodes.begin(),listNodes.end(),g_grid.m_mesh.get<Node>(p.front()))!=listNodes.end() ||
		    std::find(listNodes.begin(),listNodes.end(),g_grid.m_mesh.get<Node>(p.back()))!=listNodes.end()){
			pairNodesFace.push_back(p);
		}
	}


	std::cout<<"List pair Nodes of Face :"<<std::endl;
	for (auto n : pairNodesFace){
		for (auto i : n){
			std::cout<<i<<std::endl;
		}
	}
	if(pairNodesFace.empty()){
		std::cout<<"VIDE"<<std::endl;
	}

	return pairNodesFace;
}
std::vector<Node> Tools::getBoundaryEdge(Node firstNodeId, Node secondNodeId)
{
	auto listAllPair = Tools::getAllNodesChain(firstNodeId.id(),secondNodeId.id());
	std::vector<Node> boundaryEdge ={};
	for (auto i : listAllPair) {
		auto faceCommon=Tools::getIdOneCommonFace(i.front(),i.back());
		auto facesCommon = Tools::getFacesCommon(i.front(),i.back());
		if(facesCommon.front().id()==faceCommon &&facesCommon.back().id()==faceCommon){
			boundaryEdge.push_back(g_grid.m_mesh.get<Node>(i.front()));
			boundaryEdge.push_back(g_grid.m_mesh.get<Node>(i.back()));
			std::cout<<"Edge Boundary :"<<std::endl;
			for (auto b : boundaryEdge) {
				std::cout<<b<<std::endl;
			}
			return boundaryEdge;
		}
	}
}


void Tools::joinFaceToNodes(Node i1, Node i2)
{
	auto boundaryEdge = Tools::getBoundaryEdge(i1,i2);
	auto listFaces = Tools::getAllFacesChain(boundaryEdge.front().id(),boundaryEdge.back().id());
	auto listPairNodes = Tools::getAllNodesChain(boundaryEdge.front().id(),boundaryEdge.back().id());
	auto listFirstNodes = Tools::getListFirstNodesChain(boundaryEdge.front().id(),boundaryEdge.back().id());
	auto listSecondNodes = Tools::getListSecondNodesChain(boundaryEdge.front().id(),boundaryEdge.back().id());

	Node oldNewNodePair1;
	Node oldNewNodePair2;


	for(auto f : listFaces){
		auto pairNodesFace=Tools::getPairNodesFace(i1.id(),i2.id(),f.id());
		std::vector<int> firstPair = pairNodesFace.front();
		std::vector<int> secondPair = pairNodesFace.back();
		std::cout<< "La face : " << f.id() <<std::endl;
		std::cout<< "Pair 1 : "<<firstPair.front()<<" "<<firstPair.back()<<std::endl;
		std::cout<< "Pair 2 : "<<secondPair.front()<<" "<<secondPair.back()<<std::endl;
		auto newNodePair1= Tools::createMiddleNode(g_grid.m_mesh.get<Node>(firstPair.front()),g_grid.m_mesh.get<Node>(firstPair.back()));
		auto newNodePair2 = Tools::createMiddleNode(g_grid.m_mesh.get<Node>(secondPair.front()),g_grid.m_mesh.get<Node>(secondPair.back()));
		std::cout<< "Pair 1 new node: "<<newNodePair1<<std::endl;
		std::cout<< "Pair 2 new node: "<<newNodePair2<<std::endl;

		Node nodeFirstPair1;
		Node nodeSecondPair1;

		Node nodeFirstPair2;
		Node nodeSecondPair2;


		for (auto l : listFirstNodes){
			for (auto i : firstPair){
				if (i==l){
					nodeSecondPair1= g_grid.m_mesh.get<Node>(i);
				}
			}
			for( auto i : secondPair){
				if (i==l){
					nodeSecondPair2= g_grid.m_mesh.get<Node>(i);
				}
			}
		}
		for (auto l : listSecondNodes) {
			for (auto i : firstPair) {
				if (i == l) {
					nodeFirstPair1 = g_grid.m_mesh.get<Node>(i);
				}
			}
			for (auto i : secondPair) {
				if (i == l) {
					nodeFirstPair2 = g_grid.m_mesh.get<Node>(i);
				}
			}
		}
		std::cout<< "LES FACES : "<<std::endl;

		if(newNodePair1 == oldNewNodePair1){
			Face newFaceLeft = g_grid.m_mesh.newQuad(nodeFirstPair1,oldNewNodePair1,newNodePair2,nodeFirstPair2);
			g_grid.m_mesh.get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
			g_grid.m_mesh.get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
			g_grid.m_mesh.get<Node>(oldNewNodePair1.id()).add(newFaceLeft);
			g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceLeft);


			Face newFaceRight = g_grid.m_mesh.newQuad(nodeSecondPair1,oldNewNodePair1,newNodePair2,nodeSecondPair2);
			g_grid.m_mesh.get<Node>(nodeSecondPair1.id()).add(newFaceRight);
			g_grid.m_mesh.get<Node>(nodeSecondPair2.id()).add(newFaceRight);
			g_grid.m_mesh.get<Node>(oldNewNodePair1.id()).add(newFaceRight);
			g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceRight);
		}
		else{
			if(newNodePair1==oldNewNodePair2){
				Face newFaceLeft = g_grid.m_mesh.newQuad(nodeFirstPair1,oldNewNodePair2,newNodePair2,nodeFirstPair2);
				g_grid.m_mesh.get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
				g_grid.m_mesh.get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
				g_grid.m_mesh.get<Node>(oldNewNodePair2.id()).add(newFaceLeft);
				g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceLeft);


				Face newFaceRight = g_grid.m_mesh.newQuad(nodeSecondPair1,oldNewNodePair2,newNodePair2,nodeSecondPair2);
				g_grid.m_mesh.get<Node>(nodeSecondPair1.id()).add(newFaceRight);
				g_grid.m_mesh.get<Node>(nodeSecondPair2.id()).add(newFaceRight);
				g_grid.m_mesh.get<Node>(oldNewNodePair2.id()).add(newFaceRight);
				g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceRight);
			}
			else{
				if(newNodePair2==oldNewNodePair1){
					Face newFaceLeft = g_grid.m_mesh.newQuad(nodeFirstPair1,newNodePair1,oldNewNodePair1,nodeFirstPair2);
					g_grid.m_mesh.get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
					g_grid.m_mesh.get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
					g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceLeft);
					g_grid.m_mesh.get<Node>(oldNewNodePair1.id()).add(newFaceLeft);


					Face newFaceRight = g_grid.m_mesh.newQuad(nodeSecondPair1,newNodePair1,oldNewNodePair1,nodeSecondPair2);
					g_grid.m_mesh.get<Node>(nodeSecondPair1.id()).add(newFaceRight);
					g_grid.m_mesh.get<Node>(nodeSecondPair2.id()).add(newFaceRight);
					g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceRight);
					g_grid.m_mesh.get<Node>(oldNewNodePair1.id()).add(newFaceRight);

				}
				else{
					if(newNodePair2==oldNewNodePair2){
						Face newFaceLeft = g_grid.m_mesh.newQuad(nodeFirstPair1,newNodePair1,oldNewNodePair2,nodeFirstPair2);
						g_grid.m_mesh.get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(oldNewNodePair2.id()).add(newFaceLeft);


						Face newFaceRight = g_grid.m_mesh.newQuad(nodeSecondPair1,newNodePair1,oldNewNodePair2,nodeSecondPair2);
						g_grid.m_mesh.get<Node>(nodeSecondPair1.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(nodeSecondPair2.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(oldNewNodePair2.id()).add(newFaceRight);
					}
					else{
						Face newFaceLeft = g_grid.m_mesh.newQuad(nodeFirstPair1,newNodePair1,newNodePair2,nodeFirstPair2);
						g_grid.m_mesh.get<Node>(nodeFirstPair1.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(nodeFirstPair2.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceLeft);
						g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceLeft);


						Face newFaceRight = g_grid.m_mesh.newQuad(nodeSecondPair1,newNodePair1,newNodePair2,nodeSecondPair2);
						g_grid.m_mesh.get<Node>(nodeSecondPair1.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(nodeSecondPair2.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(newNodePair1.id()).add(newFaceRight);
						g_grid.m_mesh.get<Node>(newNodePair2.id()).add(newFaceRight);

					}
				}
			}
		}


	}
	for (auto lFace : listFaces){
		g_grid.m_mesh.deleteFace(lFace.id());
	}
}