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

        void remove(int ASheetID);

        void generateDual(GraphicView *AView);
        void generateDual(TextView *AView);
        void generateMesh(GraphicView *AView);

        void showAxis(GraphicView *AView, vtkIdType ACellID);

        void generateBlocks(GraphicView *AView);

        void surfaceMode(GraphicView *AView);

        void resetDual(GraphicView *AView);
        void resetBlock(GraphicView *AView);
        void resetBlockSmoothing(GraphicView *AView);
        void smoothBlocks(GraphicView *AView);

        void singularityGraphToggle(GraphicView *AView);
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATOREXPLORERPOLICY_H
/*----------------------------------------------------------------------------*/
