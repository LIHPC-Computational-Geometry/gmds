#include <gmds/paul/Actions_Agent.h>
using namespace gmds;

Actions::Actions(GridBuilderAround *AGrid)
   :g_grid(*AGrid), tool(&g_grid){;}

/*
Actions::~Actions(){
	delete &g_grid;
	delete this;
}*/

void Actions::deleteActions()
{
	delete this;
}

void Actions::executeDeleteFace(int faceID)
{
	g_grid.flipActivate(faceID);
}

void Actions::executeCutEdge(Node firstNodeID, Node secondNodeID)
{
	if(tool.checkExistEdge(firstNodeID.id(),secondNodeID.id(),tool.getIdOneCommonFace(firstNodeID.id(),secondNodeID.id()))){
		tool.createEdge(firstNodeID,secondNodeID);
	}
	else{
		std::cout<<"Error, no edge between nodes"<<std::endl;
	}
}

void Actions::executeGlideNode(Node node, Mesh *AMesh)
{
	auto boundaryNodes = tool.getBoundaryNodes(AMesh);
	double range;
	bool first = true;
	double newX;
	double newY;
	double newZ;
	auto listFacesNode = tool.getListFacesOfNode(node.id());
	double oldXNode=node.X();
	double oldYNode=node.Y();
	double oldZNode=node.Z();
	for (auto n : boundaryNodes){
		double newRange = tool.calcRangePoints(node,AMesh->get<Node>(n.first));
		if (first == true){
			range = newRange;
			first=false;
			newX=AMesh->get<Node>(n.first).X();
			newY=AMesh->get<Node>(n.first).Y();
			newZ=AMesh->get<Node>(n.first).Z();
		}
		else if (newRange < range){
			range = newRange;
			newX=AMesh->get<Node>(n.first).X();
			newY=AMesh->get<Node>(n.first).Y();
			newZ=AMesh->get<Node>(n.first).Z();
		}
	}

	node.X()=newX;
	node.Y()=newY;
	node.Z()=newZ;

	if(! g_grid.isValid()){
		std::cout<<"Bad Cell detected"<<std::endl;
		while(!g_grid.isValid()){
			std::cout<<"Bad Cell detected in while"<<std::endl;
			node.X()= (node.X() + oldXNode)/2;
			node.Y() = (node.Y() + oldYNode)/2;
			node.Z() = (node.Z() + oldZNode)/2;
		}
	}


}
void Actions::executeCutFace(Face AFace, int direction)
{
	if (direction==1){
		auto edgeHorizontal = tool.getHorizontalEdge(AFace);
		Actions::executeCutEdge(edgeHorizontal.front(),edgeHorizontal.back());
	}
	else if (direction == 0){
		auto edgeVertical = tool.getVerticalEdge(AFace);
		Actions::executeCutEdge(edgeVertical.front(),edgeVertical.back());
	}
	else{std::cout<<"Error input"<<std::endl;}
}

void Actions::executeGlideMaxNodeFace(Face AFace)
{
	Node ANode = tool.selectNodeMaxRange(AFace);
	Actions::executeGlideNode(ANode,g_grid.meshTarget);
}

void Actions::executeGlideMinNodeFace(Face AFace)
{
	Node ANode = tool.selectNodeMinRange(AFace);
	Actions::executeGlideNode(ANode,g_grid.meshTarget);
}
