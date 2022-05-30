#include <gmds/baptiste/action.h>

using namespace gmds;

Action::Action() {}

/*
void Action::executeAction(RLBlockSet &blockSet, int faceID)
{
	std::cout << "Entering base executeAction function ...";
}
 */

ActionDelete::ActionDelete() {}

void ActionDelete::executeAction(RLBlockSet &blockSet, int faceID)
{
	blockSet.deleteBlock(faceID);
}

ActionEdit::ActionEdit(int V, int Axis, int Range)
{
	v = !!V;
	if (Axis == 0)
	{
		axis = "x";
	}
	else
	{
		axis = "y";
	}
	range = Range;
}

void ActionEdit::executeAction(RLBlockSet &blockSet, int faceID)
{
	blockSet.editCorner(faceID, v, axis, range);
}