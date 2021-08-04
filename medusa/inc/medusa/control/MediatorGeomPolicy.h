/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIATORGEOMPOLICY_H
#define GMDS_MEDIATORGEOMPOLICY_H
/*----------------------------------------------------------------------------*/
#include <medusa/control/MediatorPolicyItf.h>
/*----------------------------------------------------------------------------*/

namespace medusa {
    class MediatorGeomPolicy: public MediatorPolicyItf {

    public:
        MediatorGeomPolicy();
        virtual ~MediatorGeomPolicy();

        /*
         * Pick a geom surface
         */
        int selectVTKCell(GraphicView* AView,
                          const int ADim, const vtkIdType AID, int AAxis);
        int selectVTKCell(TextView* AView,
                          const int ADim, const vtkIdType AID);
        void help(TextView* AView);

        void remove(GraphicView *AView);

        /*
         * Will generate the frame field if there is no one, will generate the geom cut otherwise.
         */
        void generate(GraphicView *AView);
        void generate(TextView *AView);

        /*
         * Remove the last geom surface selected
         */
        void removeItem(int ASheetID, GraphicView *AView);
        void removeItem(int ASheetID, TextView *AView);

        //----------------------------------------------
        /*
         * Do nothing for now
         */
        void showAxis(GraphicView *AView, vtkIdType ACellID);

        void surfaceMode(GraphicView *AView);
        //----------------------------------------------

        /*
         * Undo the last "generate" action done
         */
        void undoGenerate(GraphicView *AView);

        void singularityGraphToggle(GraphicView *AView);

        /*
         * Calling this method switch the opacity of the cutted zone alternatively between transparent to hided to opaque
         */
        void opacity(GraphicView *AView);

        /*
         * toggle the geom surfaces colors
         */
        void color(GraphicView *AView);

    protected:

        /*
         * 0 = input mode : Frame field not generated
         * 1 = frame mode : Frame field generated, can do all frame operations
         * 2 = geom mode  : can do all geom operation
         */
        int m_mode = 0;

    };

}


#endif //GMDS_MEDIATORGEOMPOLICY_H
