//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMAND_H
#define GMDS_MEDUSACOMMAND_H

#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include "MediatorPolicyItf.h"

namespace medusa{
    class MedusaCommand{
    public:
        enum State {SUCCES,FAILURE};
    public:
        MedusaCommand(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView){
            m_policy = APolicy;
            m_graphic_collaborator = AGraphicView;
            m_text_collaborator = ATextView;
        }
        virtual ~MedusaCommand(){};
        virtual State execute() = 0;
        virtual State undo() = 0;

    protected:
        MediatorPolicyItf* m_policy;
        GraphicView* m_graphic_collaborator;
        TextView* m_text_collaborator;
    };
}

#endif //GMDS_MEDUSACOMMAND_H
