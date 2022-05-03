#include "gmds/paul/Grid.h"
#include <gmds/igalgo/GridBuilder.h>
using namespace gmds;

Grid::Grid(const int AX, const int AY)
  :m_X(AX), m_Y(AY){;}

int Grid::getX() const
{return  m_X;}

int Grid::getY() const
{return m_Y;}

GridBuilderAround::GridBuilderAround(Mesh *AMesh, const TInt ADim)
   :m_mesh(*AMesh), m_dim(ADim){;}

const Mesh
GridBuilderAround::getMesh() const
{return  m_mesh;}

int GridBuilderAround::getDim() const
{return m_dim;}

/*int GridBuilderAround::getActivate(const int faceID)
{
      //return m_mesh.get<Face>(faceID).;
}*/

GridBuilderAround::~GridBuilderAround()
{}

/*----------------------------------------------------------------------------*/
void GridBuilderAround::executeGrid2D(const gmds::TInt ANb) {
	double Xmax = 0;
	double Ymax = 0;
	double Xmin = 0;
	double Ymin = 0;
	for(auto n : m_mesh.nodes()){
		if ( m_mesh.get<Node>(n).X()>Xmax){
			Xmax=m_mesh.get<Node>(n).X();
		}
		if ( m_mesh.get<Node>(n).Y()>Ymax){
			Ymax=m_mesh.get<Node>(n).Y();
		}
		if ( m_mesh.get<Node>(n).X()<Xmin){
			Xmin=m_mesh.get<Node>(n).X();
		}
		if ( m_mesh.get<Node>(n).Y()<Ymin){
			Ymin=m_mesh.get<Node>(n).Y();
		}
	}
	std::cout<<"Valeur X max :"<<Xmax<<std::endl;
	std::cout<<"Valeur Y max :"<<Ymax<<std::endl;
	std::cout<<"Valeur X min :"<<Xmin<<std::endl;
	std::cout<<"Valeur Y min :"<<Ymin<<std::endl;
	double rangeX = abs(Xmin-Xmax);
	double rangeY = abs(Ymin-Ymax);
	std::cout<<"Valeur range X :"<<rangeX<<std::endl;
	std::cout<<"Valeur range Y :"<<rangeY<<std::endl;
	m_mesh.clear();
	gridBuild2D(ANb, rangeX, ANb, rangeY,Xmin,Ymin);
}

void GridBuilderAround::gridBuild2D(const gmds::TInt AXNb,
                     const gmds::TCoord AXStep,
                     const gmds::TInt AYNb,
                     const gmds::TCoord AYStep, const gmds::TCoord Xmin, const gmds::TCoord Ymin)
{

	std::vector<TCellID> node_ids;
	const gmds::TInt N = AXNb * AYNb;
	node_ids.reserve(N);

	for (auto x = 0; x < AXNb; x++) {
		for (auto y = 0; y < AYNb; y++) {
			Node n = m_mesh.newNode(x*AXStep/(AXNb-1), y*AYStep/(AYNb-1), 0);
			node_ids.push_back(n.id());
		}
	}
	for (auto k = 0; k < N - AYNb; k++) {
		if ((k + 1) % AYNb == 0) {
			continue;
		}

		Face q = m_mesh.newQuad(node_ids.at(k),
										node_ids.at(k + AYNb),			// [x + 1][y],
										node_ids.at(k + AYNb + 1),		// [x + 1][y + 1],
										node_ids.at(k + 1));			// [x][y + 1],
		m_mesh.get<Node>(node_ids.at(k)).add(q);
		m_mesh.get<Node>(node_ids.at(k + AYNb)).add(q);
		m_mesh.get<Node>(node_ids.at(k + AYNb + 1)).add(q);
		m_mesh.get<Node>(node_ids.at(k + 1)).add(q);
	}
	//add variable activate at faces
	for (auto face_id:m_mesh.faces()){
		Face f=m_mesh.get<Face>(face_id);
		std::vector<Node> f_nodes = f.get<Node>();
		for (auto n : f_nodes){
			activate->set(f.id(),1); //attribution val aux faces
		}
	}

}

int GridBuilderAround::getActivate(gmds::Face AFace) {
	int  valueActivateFace = m_mesh.getVariable<int,GMDS_FACE>("exist")->value(AFace.id());
	return valueActivateFace;
}


void GridBuilderAround::flipActivate(const int faceID)
{
	activate->set(faceID,0);
}
