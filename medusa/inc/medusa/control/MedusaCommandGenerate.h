//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMANDGENERATE_H
#define GMDS_MEDUSACOMMANDGENERATE_H


#include "MedusaCommand.h"

namespace medusa{
    class MedusaCommandGenerate: public MedusaCommand{

    public:
        MedusaCommandGenerate(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView);
        State execute();
        State undo();

    protected:

    };
}


#endif //GMDS_MEDUSACOMMANDGENERATE_H
