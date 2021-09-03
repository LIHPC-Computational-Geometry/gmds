/*----------------------------------------------------------------------------*/
#include "medusa/control/MedusaCommandBlockSmoothing.h"
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
MedusaCommandBlockSmoothing::MedusaCommandBlockSmoothing(MediatorPolicyItf* APolicy, GraphicView* AGraphicView, TextView* ATextView): MedusaCommand(APolicy, AGraphicView, ATextView) {

}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandBlockSmoothing::execute() {
    m_policy->smoothBlocks(m_graphic_collaborator);
    return SUCCES;
}
/*----------------------------------------------------------------------------*/
MedusaCommand::State MedusaCommandBlockSmoothing::undo() {
    m_policy->resetBlockSmoothing(m_graphic_collaborator);
}
/*----------------------------------------------------------------------------*/
