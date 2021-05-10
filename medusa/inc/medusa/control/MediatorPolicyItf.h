/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIATORPOLICYITF_H
#define GMDS_MEDIATORPOLICYITF_H
/*----------------------------------------------------------------------------*/
#include <vtkType.h>
/*----------------------------------------------------------------------------*/
namespace medusa {
    /*------------------------------------------------------------------------*/
    class GraphicView;
    class TextView;
    /*------------------------------------------------------------------------*/
    class MediatorPolicyItf {
    public:
        virtual int selectVTKCell(GraphicView* AView,
                                   const int ADim, const vtkIdType AID, int AAxis)=0;
        virtual int selectVTKCell(TextView* AView,
                                  const int ADim, const vtkIdType AID)=0;
        virtual void help(TextView* AView) = 0;
        virtual void remove(int ASheetID) = 0;
        virtual void generateDual(GraphicView *AView) = 0;
        virtual void generateDual(TextView *AView) = 0;

        virtual void showAxis(GraphicView *AView, vtkIdType ACellID) = 0;

        virtual void generateBlocks(GraphicView *AView) = 0;

        virtual void surfaceMode(GraphicView *AView) = 0;

        virtual void resetDual(GraphicView *AView) = 0;

        virtual void resetBlock(GraphicView *AView) = 0;

        virtual void smoothBlocks(GraphicView *AView) = 0;
        virtual void resetBlockSmoothing(GraphicView *AView) = 0;

        virtual void singularityGraphToggle(GraphicView *AView) = 0;
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATORPOLICYITF_H
/*----------------------------------------------------------------------------*/
