//
// Created by calderans on 21/07/22.
//

#include <gmds/blockMesher/BlockClassificator.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BlockClassificator::BlockClassificator(gmds::Mesh *ABlocks, gmds::cad::GeomMeshLinker *ALinker, gmds::cad::FACManager *AManager):
   m_blocks(ABlocks),
   m_linker(ALinker),
   m_manager(AManager) {

}
/*----------------------------------------------------------------------------*/
BlockClassificator::~BlockClassificator() {}
/*----------------------------------------------------------------------------*/
void BlockClassificator::blockCreation(int AMTrix[5][5][5]){
	int S[6][6][6];

	for(int k = 0; k<5; k++){
		for(int j = 0; j<5; j++){
			for(int i = 0; i<5; i++){
				std::cout<<AMTrix[i][j][k]<<" ";
				if(AMTrix[i][j][k] == 1){
					S[i][j][k] = 1;
					S[i+1][j][k] = 1;
					S[i+1][j+1][k] = 1;
					S[i][j+1][k] = 1;
					S[i][j][k+1] = 1;
					S[i+1][j][k+1] = 1;
					S[i+1][j+1][k+1] = 1;
					S[i][j+1][k+1] = 1;
				}
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}

	int geom_scale = 5;

	for(int k = 0; k<6; k++){
		for(int j = 0; j<6; j++){
			for(int i = 0; i<6; i++){
				if(S[i][j][k] == 1){
					S[i][j][k] = 1;
					Node n = m_blocks->newNode(i*geom_scale,j*geom_scale,k*geom_scale);
					S[i][j][k] = n.id();
				}else{
					S[i][j][k] = -1;
				}
			}
		}
	}

	for(int k = 0; k<5; k++){
		for(int j = 0; j<5; j++){
			for(int i = 0; i<5; i++){
				if(AMTrix[i][j][k] == 1){
					m_blocks->newHex(S[i][j][k],S[i+1][j][k],S[i+1][j+1][k],S[i][j+1][k],
					                 S[i][j][k+1],S[i+1][j][k+1],S[i+1][j+1][k+1],S[i][j+1][k+1]);
				}
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
void BlockClassificator::blockClassification(){

	std::cout<<"Auto classification"<<std::endl;

	for(auto n : m_blocks->nodes()){

		Node node = m_blocks->get<Node>(n);
		math::Point n_p = node.point();

		int dim,id;

		getEntityToClassify(n_p,dim,id);

		classifyBlockEntity(n, 0, dim, id);
	}

	for(auto e : m_blocks->edges()){
		Edge edge = m_blocks->get<Edge>(e);

		int nbBndFaces = 0;
		for(auto f : edge.get<Face>()){
			if(f.get<Region>().size() == 1)
				nbBndFaces++;
		}
		if (nbBndFaces == 2){
			math::Point center = edge.center();

			int dim,id;

			getEntityToClassify(center,dim,id);

			classifyBlockEntity(e, 1, dim, id);
		}
	}

	for(auto f : m_blocks->faces()){
		Face face = m_blocks->get<Face>(f);

		if(face.get<Region>().size() == 1){

			math::Point center = face.center();

			int dim,id;

			getEntityToClassify(center,dim,id);

			classifyBlockEntity(f, 2, dim, id);
		}
	}

}
/*----------------------------------------------------------------------------*/
void BlockClassificator::getEntityToClassify(math::Point APoint, int &ADim, int &AID){
	int dim = -1;
	int id = -1;
	double distance = 9999999;
	std::vector<cad::GeomSurface*> surfaces;
	std::vector<cad::GeomCurve*> curves;
	std::vector<cad::GeomPoint*> points;

	m_manager->getSurfaces(surfaces);
	for (auto s : surfaces) {
		math::Point project(APoint.X(),APoint.Y(),APoint.Z());
		s->project(project);
		if(APoint.distance(project) <= distance){
			distance = APoint.distance(project);
			dim = 2;
			id = s->id();
		}
	}
	m_manager->getCurves(curves);
	for (auto c : curves) {
		math::Point project(APoint.X(),APoint.Y(),APoint.Z());
		c->project(project);
		if(APoint.distance(project) <= distance){
			distance = APoint.distance(project);
			dim = 1;
			id = c->id();
		}
	}
	m_manager->getPoints(points);
	for (auto p : points) {
		math::Point project(APoint.X(),APoint.Y(),APoint.Z());
		p->project(project);
		if(APoint.distance(project) <= distance){
			distance = APoint.distance(project);
			dim = 0;
			id = p->id();
		}
	}
	ADim = dim;
	AID = id;
}
/*----------------------------------------------------------------------------*/
void BlockClassificator::classifyBlockEntity(TCellID ACellID, int ACellDim, int ADim, int AID){

	Node node;
	math::Point project;
	Edge edge;
	Face face;

	if(ACellDim == 0) {
		node = m_blocks->get<Node>(ACellID);
		project.setXYZ(node.point().X(), node.point().Y(), node.point().Z());
	}else if(ACellDim == 1){
		edge = m_blocks->get<Edge>(ACellID);
	}else if(ACellDim == 2){
		face = m_blocks->get<Face>(ACellID);
	}

	if(ADim == 0){
		m_linker->linkToPoint(node,AID);
		cad::GeomPoint* geomPoint = m_manager->getPoint(AID);
		geomPoint->project(project);
	}else if(ADim == 1){
		if(ACellDim == 0){
			m_linker->linkToCurve(node,AID);
			cad::GeomCurve* geomCurve = m_manager->getCurve(AID);
			geomCurve->project(project);
		}else if(ACellDim == 1){
			m_linker->linkToCurve(edge,AID);
		}
	}else if(ADim == 2){
		if(ACellDim == 0){
			m_linker->linkToSurface(node, AID);
			cad::GeomSurface* geomSurface = m_manager->getSurface(AID);
			geomSurface->project(project);
		}else if(ACellDim == 1){
			m_linker->linkToSurface(edge, AID);
		}else if(ACellDim == 2){
			m_linker->linkToSurface(face, AID);
		}

	}

	if(ACellDim == 0){
		node.setPoint(project);
	}
}
/*----------------------------------------------------------------------------*/
Mesh* BlockClassificator::getBlocks(){
	return m_blocks;
}
