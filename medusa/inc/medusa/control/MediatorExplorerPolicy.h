/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIATOREXPLORERPOLICY_H
#define GMDS_MEDIATOREXPLORERPOLICY_H
/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorPolicyItf.h>
/*----------------------------------------------------------------------------*/
namespace medusa {
    class MediatorExplorerPolicy : public MediatorPolicyItf {

    public:
        MediatorExplorerPolicy();
        virtual ~MediatorExplorerPolicy();
        int selectVTKCell(GraphicView* AView,
                           const int ADim, const vtkIdType AID, int AAxis);
        int selectVTKCell(TextView* AView,
                          const int ADim, const vtkIdType AID);
        void help(TextView* AView);

        void remove(GraphicView *AView);

        void removeItem(int ASheetID, GraphicView *AView);
        void removeItem(int ASheetID, TextView *AView);

        void generate(GraphicView *AView);
        void generate(TextView *AView);

        void showAxis(GraphicView *AView, vtkIdType ACellID);

        void surfaceMode(GraphicView *AView);

        void undoGenerate(GraphicView *AView);

        void singularityGraphToggle(GraphicView *AView);

        void opacity(GraphicView *AView);

        void color(GraphicView *AView);
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATOREXPLORERPOLICY_H
/*----------------------------------------------------------------------------*/
