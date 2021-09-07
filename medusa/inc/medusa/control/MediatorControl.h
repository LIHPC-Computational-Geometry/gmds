/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIATORCONTROL_H
#define GMDS_MEDIATORCONTROL_H
/*----------------------------------------------------------------------------*/
#include <medusa/view/MedusaObserver.h>
#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include <medusa/control/MediatorPolicyItf.h>
/*----------------------------------------------------------------------------*/
#include <vtkType.h>
#include "MedusaCommand.h"
#include "MedusaCommandPick.h"
/*----------------------------------------------------------------------------*/
namespace medusa{
    class MediatorControl:public MedusaObserver{
    public:

        static MediatorControl* getInstance();
        static void initialize();
        void update();

        void addCollaborator(GraphicView* AView);
        void addCollaborator(TextView* AView);

        void setPolicy(MediatorPolicyItf* APolicy);

        void select(int ADIm, vtkIdType ACellID);
        void help();

        void generate();
        void toggleDual();

        void wireframe();

        void opacity();

        void remove();

        int removeActor(vtkSmartPointer<vtkActor> ADeletedActor);

        void removeItem(int AItemID);

        void surfaceMode();
        void finalMesh();
        void undo();

        void color();

        void singularityGraphToggle();

        void changePolicy();

    protected:

        MediatorControl();
        virtual ~MediatorControl();

    private:
        static MediatorControl* m_instance;

        MediatorPolicyItf* m_policy;
        std::vector<GraphicView*>   m_graphic_collaborator;
        std::vector<TextView*>      m_text_collaborator;
        std::vector<MedusaCommand*> m_commands;

        //flag for the visualisation of the dual, -1 the dual isnt generated yet, 0 the dual is hidden and 1 the dual is visible
        int m_dual_state = -1;

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATORCONTROL_H
/*----------------------------------------------------------------------------*/
