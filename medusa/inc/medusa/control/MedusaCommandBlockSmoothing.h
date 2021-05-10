//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMAND_MESH_H
#define GMDS_MEDUSACOMMAND_MESH_H


#include "MedusaCommand.h"

namespace medusa{
    class MedusaCommandBlockSmoothing: public MedusaCommand{

    public:
        MedusaCommandBlockSmoothing(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView);
        State execute();
        State undo();

    protected:

    };
}


#endif //GMDS_MEDUSACOMMANDBLOCK_H
