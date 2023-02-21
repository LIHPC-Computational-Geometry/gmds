//
// Created by bourmaudp on 02/12/22.
//
#include <gmds/blocking_Quality/BlockingQuality.h>

namespace gmds{

double blockQuality(Region ARegion)
{
	double bQuality=0;
	std::vector<Node> listNode =ARegion.getAll<Node>();

	quality::HexQuality he = quality::HexQuality::build(listNode[0].point(),listNode[1].point(),
	                                                    listNode[2].point(),listNode[3].point(),
	                                                    listNode[4].point(),listNode[5].point(),
	                                                    listNode[6].point(),listNode[7].point());


	bQuality= (he.diagonal() + he.edgeRatioModified() + he.scaledJacobian())/3;
	return bQuality;
}

double allBlocksQuality(Mesh *AMesh){
	double gBlocksQuality = 0;
	int nbRegions = AMesh->getNbRegions();

	for (auto r : AMesh->regions()){
		gBlocksQuality = gBlocksQuality+ blockQuality(AMesh->get<Region>(r));
	}

	gBlocksQuality = gBlocksQuality/nbRegions;
	return gBlocksQuality;
}

double layeringNodes(Mesh *AMesh, const Mesh *AImprintMesh){
	int nbNodesCommon = 0;
	for(auto nAMesh : AMesh->nodes()){
		for (auto nAImprintMesh : AImprintMesh->nodes()){
			if(AMesh->get<Node>(nAMesh).X() == AImprintMesh->get<Node>(nAImprintMesh).X()  &&
			    AMesh->get<Node>(nAMesh).Y() == AImprintMesh->get<Node>(nAImprintMesh).Y() &&
			    AMesh->get<Node>(nAMesh).Z() == AImprintMesh->get<Node>(nAImprintMesh).Z() ){
				nbNodesCommon++;
			}
		}
	}

	double featRatio = (double)nbNodesCommon/(double)AImprintMesh->getNbNodes();
	//std::cout<<"FEAT RATIO "<<featRatio<<std::endl;
	return featRatio;


}

double boundaryEdgeValence(Mesh *AMesh){
	double valenceBoundary = 0;
	for(auto e : AMesh->edges()) {
		if (AMesh->getVariable<int, GMDS_EDGE>("boundaryEdge")->value(e) == 1) {
			if (AMesh->get<Edge>(e).nbRegions() == 1) {
				valenceBoundary += 0.25;
			}
			else if (AMesh->get<Edge>(e).nbRegions() == 3) {
				valenceBoundary -= 0.25;
			}
			else if (AMesh->get<Edge>(e).nbRegions() == 4) {
				valenceBoundary -= 0.5;
			}
		}
	}
	//std::cout<<"Valence Region "<<valenceBoundary<<std::endl;
	return valenceBoundary;
}
double internEdgeValence(Mesh *AMesh){
	double valenceIntern = 0;
	for(auto e : AMesh->edges()) {
		if (AMesh->getVariable<int, GMDS_EDGE>("boundaryEdge")->value(e) == 0) {
			if (AMesh->get<Edge>(e).nbRegions() == 3) {
				valenceIntern += 0.25;
			}
			else if (AMesh->get<Edge>(e).nbRegions() == 5) {
				valenceIntern -= 0.25;
			}
		}
	}
	//std::cout<<"Valence Intern "<<valenceIntern<<std::endl;
	return valenceIntern;
}

double globalEdgeValence(Mesh *AMesh){
	double valenceGlobal = 0;
	valenceGlobal = internEdgeValence(AMesh) + boundaryEdgeValence(AMesh);
	//std::cout<<"Valence global "<<valenceGlobal<<std::endl;
	return valenceGlobal;
}

double ratioGlobalEdgeValence(Mesh *AMesh){
	return globalEdgeValence(AMesh)/AMesh->getNbEdges();
}

int blockingSimplicity(Mesh *AMesh){
	return AMesh->getNbRegions();
}

double blockingQuality(Mesh *AMesh, Mesh *AGeometry){
	double blockingQuality = 0;
	double blockingQualityWithoutGeo = 0;

	blockingQualityWithoutGeo = (allBlocksQuality(AMesh) + ratioGlobalEdgeValence(AMesh)) / double(blockingSimplicity(AMesh));

	blockingQuality = (blockingQualityWithoutGeo + layeringNodes(AMesh,AGeometry));

	return blockingQuality;
}

std::unordered_map<TCellID ,std::vector<TCellID>> mapNodes2Points(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapN2P;
	for(auto n : ABlocks->nodes()){
		TCellID p = ALinker->getGeomId(ABlocks->get<Node>(n));
		mapN2P[p].push_back(n);
	}
/*
	for (auto e : mapN2P){
		std::cout<<"Le Point : "<<e.first<<" avec les noeuds :"<<std::endl;
		for (auto n : e.second){
			std::cout<<n<<";";
		}
		std::cout<<""<<std::endl;
	}*/
	return mapN2P;
}


std::unordered_map<TCellID ,std::vector<TCellID>> mapEdges2Curves(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapN2P;
	for(auto e : ABlocks->edges()){
		TCellID c = ALinker->getGeomId(ABlocks->get<Edge>(e));
		mapN2P[c].push_back(e);
	}

	for (auto e : mapN2P){
		std::cout<<"Le Courbes : "<<e.first<<" avec les noeuds :"<<std::endl;
		for (auto n : e.second){
			std::cout<<n<<";";
		}
		std::cout<<""<<std::endl;
	}
	return mapN2P;
}

std::unordered_map<TCellID ,std::vector<TCellID>> mapFaces2Surfaces(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapN2P;
	for(auto f : ABlocks->faces()){
		TCellID s = ALinker->getGeomId(ABlocks->get<Face>(f));
		mapN2P[s].push_back(f);
	}

	for (auto e : mapN2P){
		std::cout<<"Le Surfaces : "<<e.first<<" avec les infos :"<<std::endl;
		for (auto n : e.second){
			std::cout<<n<<";";
		}
		std::cout<<""<<std::endl;
	}
	return mapN2P;
}

bool linkNodes2Points(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapN2P= mapNodes2Points(ABlocks,AGeom,ALinker);
	for (auto n : AGeom->nodes()){
		if( mapN2P[n].empty()){
			return false;
		}
	}
	return true;



}

bool linkEdges2Curves(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapE2C= mapEdges2Curves(ABlocks,AGeom,ALinker);
	for (auto c : AGeom->edges()){
		if( mapE2C[c].empty()){
			return false;
		}
	}
	return true;
}

bool linkFaces2Surfaces(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapF2S= mapFaces2Surfaces(ABlocks,AGeom,ALinker);
	for (auto s : AGeom->faces()){
		if( mapF2S[s].empty()){
			return false;
		}
	}
	return true;
}


std::vector<Node> noLinkNodes2Points(Mesh *ABlocks, Mesh *AGeom, cad::GeomMeshLinker *ALinker){
	std::unordered_map<TCellID ,std::vector<TCellID>> mapN2P= mapNodes2Points(ABlocks,AGeom,ALinker);
	std::vector<Node> noLinkPoints;
	for (auto n : AGeom->nodes()){
		if( mapN2P[n].empty()){
			noLinkPoints.push_back(AGeom->get<Node>(n));
		}
	}
	std::cout<<"List noeuds non link : "<<std::endl;
	for(auto n : noLinkPoints){
		std::cout<<"Noeud : "<<n<<std::endl;
	}
	return noLinkPoints;
}


}

