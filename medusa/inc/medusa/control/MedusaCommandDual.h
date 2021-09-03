//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMANDDUAL_H
#define GMDS_MEDUSACOMMANDDUAL_H


#include "MedusaCommand.h"

namespace medusa{
    class MedusaCommandDual: public MedusaCommand{

    public:
        MedusaCommandDual(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView);
        State execute();
        State undo();

    protected:

    };
}


#endif //GMDS_MEDUSACOMMANDDUAL_H
