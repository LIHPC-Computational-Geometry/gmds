#include <gmds/paul/Actions_Agent.h>
using namespace gmds;

Actions::Actions(GridBuilderAround *AGrid)
   :g_grid(*AGrid){;}

void Actions::executeDeleteFace(const int faceID)
{
	g_grid.flipActivate(faceID);
}

void Actions::executeCutEdge(gmds::Edge edgeID)
{

}

void Actions::executeGlideNode(gmds::Node nodeID)
{

}
