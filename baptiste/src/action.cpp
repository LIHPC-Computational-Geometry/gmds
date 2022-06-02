#include <gmds/baptiste/action.h>

using namespace gmds;

Action::Action() {}

/*
void Action::executeAction(RLBlockSet &blockSet, int faceID)
{
	std::cout << "Base executeAction" << "\n";

	if (type == "delete")
	{
		ActionDelete* child = dynamic_cast<ActionDelete*>(this);
		child->executeAction(blockSet, faceID);
	}
	else if (type == "edit")
	{
		ActionEdit* child = dynamic_cast<ActionEdit*>(this);
		child->executeAction(blockSet, faceID);
	}

	if (dynamic_cast<ActionDelete*>(this) == nullptr)
	{
		ActionEdit* child = dynamic_cast<ActionEdit*>(this);
		child->executeAction(blockSet, faceID);
	}
	else if (dynamic_cast<ActionEdit*>(this) == nullptr)
	{
		ActionDelete* child = dynamic_cast<ActionDelete*>(this);
		child->executeAction(blockSet, faceID);
	}
}
 */

ActionDelete::ActionDelete() {}

void ActionDelete::executeAction(RLBlockSet &blockSet, int faceID)
{
	std::cout << "Executing action of type delete" << "\n";
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
	std::cout << "Executing action of type edit" << "\n";
	blockSet.editCorner(faceID, v, axis, range);
	//blockSet.editCornerBis(faceID, v, axis, range);
}