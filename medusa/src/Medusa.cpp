/*----------------------------------------------------------------------------*/
#include <medusa/Medusa.h>
#include <medusa/model/MedusaBackEnd.h>
#include <medusa/view/GraphicView.h>
#include <medusa/control/MediatorExplorerPolicy.h>
#include <medusa/control/MediatorDuBloPolicy.h>
/*----------------------------------------------------------------------------*/
#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkHexahedron.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
using namespace gmds;
/*----------------------------------------------------------------------------*/
Medusa::Medusa(const Medusa::Mode AMode) {
    init(AMode);
}
/*----------------------------------------------------------------------------*/
Medusa::~Medusa() {
    for(auto v: m_gviews){
        if(v!=NULL)
            delete v;
    }
    for(auto v: m_tviews){
        if(v!=NULL)
            delete v;
    }
}
/*----------------------------------------------------------------------------*/
void Medusa::setMode(const Medusa::Mode AMode) {

    MediatorPolicyItf* policy=NULL;
    if(AMode==Medusa::DUAL_BLOCKING)
        policy = new MediatorDubloPolicy();
    else if(AMode==Medusa::EXPLORATION)
        policy  = new MediatorExplorerPolicy();
    else
        throw GMDSException("Medusa Error: unknown mode!");

    MediatorControl::getInstance()->setPolicy(policy);
}
/*----------------------------------------------------------------------------*/
void Medusa::launch() {
    if(m_last_gview==NULL)
        throw GMDSException("Medusa requires at least one graphic view to be launched!");
    m_last_gview->start();
}
/*----------------------------------------------------------------------------*/
void Medusa::load(const std::string &AFileName) {
    MedusaBackend::getInstance()->loadVTK(AFileName);
}
/*----------------------------------------------------------------------------*/
void Medusa::init(const Medusa::Mode AMode) {
    //Important, otherwise the MedusaBackend instance is not instanciated!!!
    MedusaBackend::initialize();
    MediatorControl::initialize();
    setMode(AMode);
    MedusaBackend::getInstance()->addObserver(MediatorControl::getInstance());

}
/*----------------------------------------------------------------------------*/
GraphicView* Medusa::newGraphicView(const std::string& AName,
                                    GraphicView::ViewType AType){
    GraphicView* v = new GraphicView(AName, AType);
    m_gviews.push_back(v);
    MedusaBackend::getInstance()->addObserver(v);
    MediatorControl::getInstance()->addCollaborator(v);
    m_last_gview=v;
    return v;
}
/*----------------------------------------------------------------------------*/
TextView* Medusa::newTextView() {
    TextView* v = new TextView();
    m_tviews.push_back(v);
    MedusaBackend::getInstance()->addObserver(v);
    MediatorControl::getInstance()->addCollaborator(v);
    return v;
}
/*----------------------------------------------------------------------------*/
