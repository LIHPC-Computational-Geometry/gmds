//
// Created by calderans on 28/01/20.
//

#ifndef GMDS_MEDUSACOMMANDSHEET_H
#define GMDS_MEDUSACOMMANDSHEET_H


#include "MedusaCommand.h"

namespace medusa{
    class MedusaCommandSheet: public MedusaCommand{

    public:
        MedusaCommandSheet(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView,int ADIM, vtkIdType ACellID);
        State execute();
        State undo();

    protected:
        int m_surfID;
        int m_DIM;
        vtkIdType m_CellID;
    };
}


#endif //GMDS_MEDUSACOMMANDSHEET_H
