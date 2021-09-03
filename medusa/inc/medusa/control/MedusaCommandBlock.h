//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMANDBLOCK_H
#define GMDS_MEDUSACOMMANDBLOCK_H


#include "MedusaCommand.h"

namespace medusa{
    class MedusaCommandBlock: public MedusaCommand{

    public:
        MedusaCommandBlock(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView);
        State execute();
        State undo();

    protected:

    };
}


#endif //GMDS_MEDUSACOMMANDBLOCK_H
