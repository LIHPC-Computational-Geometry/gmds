#include <gmds/baptiste/action.h>

using namespace gmds;

Action::Action() {}

ActionDelete::ActionDelete() {}

void ActionDelete::executeAction(RLBlockSet &blockSet, int faceID)
{
	blockSet.deleteBlock(faceID);
}

ActionEdit::ActionEdit(int iV, int iAxis, int iRange)
{
	v = !!iV;
	if (iAxis == 0)
	{
		axis = "x";
	}
	else
	{
		axis = "y";
	}
	range = iRange;
}

void ActionEdit::executeAction(RLBlockSet &blockSet, int faceID)
{
	blockSet.editCorner(faceID, v, axis, range);
}