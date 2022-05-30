#ifndef GMDS_ACTION_H
#define GMDS_ACTION_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BAPTISTE_export.h"
#include <gmds/ig/Mesh.h>
#include "gmds/baptiste/RLBlockSet.h"
/*----------------------------------------------------------------------------*/

namespace gmds
{
	class Action
   {
	 public:
	   Action();

	   virtual void executeAction(RLBlockSet &blockSet, int faceID) = 0;
   };

	class ActionDelete : public Action
   {
	 public:
	   ActionDelete();

	   void executeAction(RLBlockSet &blockSet, int faceID) override;
   };

   class ActionEdit : public Action
   {
	 public:
	   ActionEdit(int vInt, int axisInt, int range);

	   bool v;
	   std::string axis;
	   int range;

	   void executeAction(RLBlockSet &blockSet, int faceID) override;
   };
}
#endif     // GMDS_ACTION_H