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

        /*
         * Pick a cell to start the dual sheet creation
         */
        int selectVTKCell(GraphicView* AView,
                           const int ADim, const vtkIdType AID, int AAxis);
        int selectVTKCell(TextView* AView,
                          const int ADim, const vtkIdType AID);
        void help(TextView* AView);

        void remove(GraphicView *AView);

        /*
         * First, generate the dual structure, if the dual structure is already generated, generate the block structure instead
         */
        void generate(GraphicView *AView);
        void generate(TextView *AView);

        /*
         * Remove a dual sheet of ID AsheetID
         */
        void removeItem(int ASheetID, GraphicView *AView);
        void removeItem(int ASheetID, TextView *AView);

        /*
         * Show the 3 frame field direction at the cell of ID ACellID in the graphic view AView
         */
        void showAxis(GraphicView *AView, vtkIdType ACellID);

        /*
         * Toggle the surface mode for the dual sheet, instead of viewing tetras, only the plane used to create the sheet is shown
         */
        void surfaceMode(GraphicView *AView);

        /*
         * Undo the last generate action done
         */
        void undoGenerate(GraphicView *AView);

        void singularityGraphToggle(GraphicView *AView);

        void opacity(GraphicView *AView);

        void color(GraphicView *AView);

    protected:
        int m_mode = 0;
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDIATORDUBLOPOLICY_H
/*----------------------------------------------------------------------------*/
