/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIATORDUBLOPOLICY_H
#define GMDS_MEDIATORDUBLOPOLICY_H
/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorPolicyItf.h>
/*----------------------------------------------------------------------------*/
namespace medusa {
    class MediatorDubloPolicy: public MediatorPolicyItf {

    public:
        MediatorDubloPolicy();
        virtual ~MediatorDubloPolicy();
        int selectVTKCell(GraphicView* AView,
                           const int ADim, const vtkIdType AID, int AAxis);
        int selectVTKCell(TextView* AView,
                          const int ADim, const vtkIdType AID);
        void help(TextView* AView);

        void generateDual(GraphicView *AView);
        void generateDual(TextView *AView);

        void smoothBlocks(GraphicView *AView);

        void remove(int ASheetID);

        void showAxis(GraphicView *AView, vtkIdType ACellID);

        void generateBlocks(GraphicView *AView);

        void surfaceMode(GraphicView *AView);

        void resetDual(GraphicView *AView);

        void resetBlock(GraphicView *AView);
        void resetBlockSmoothing(GraphicView *AView);

        void singularityGraphToggle(GraphicView *AView);
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATORDUBLOPOLICY_H
/*----------------------------------------------------------------------------*/
