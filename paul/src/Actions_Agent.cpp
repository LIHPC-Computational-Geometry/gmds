#include <gmds/paul/Actions_Agent.h>
using namespace gmds;

Actions::Actions(GridBuilderAround *AGrid)
   :g_grid(*AGrid), tool(&g_grid){;}

void Actions::executeDeleteFace(const int faceID)
{
	g_grid.flipActivate(faceID);
}

void Actions::executeCutEdge(Node firstNodeID, Node secondNodeID)
{
	if(tool.checkExistEdge(firstNodeID.id(),secondNodeID.id(),tool.getIdOneCommonFace(firstNodeID.id(),secondNodeID.id()))){
		tool.joinFaceToNodes(firstNodeID,secondNodeID);
	}
	else{
		std::cout<<"Error, no edge between nodes"<<std::endl;
	}
}

void Actions::executeGlideNode(gmds::Node nodeID)
{

}
