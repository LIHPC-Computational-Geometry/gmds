/*----------------------------------------------------------------------------*/
#ifndef GMDS_TEXTVIEW_H
#define GMDS_TEXTVIEW_H
/*----------------------------------------------------------------------------*/
#include <medusa/view/View.h>
/*----------------------------------------------------------------------------*/
namespace medusa{
    class TextView: public View {
    public:
        TextView();
        void selectVTKCell(int ADIm, vtkIdType ACellID);
        void refresh();
        void update();
        void remove(int ASheetID);

        void createSurface(const std::vector<vtkIdType> AIDs, int ASheetID);
        void help();
        void dualError();

        void generateDual();
        int getInputAxis();
    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_TEXTVIEW_H
/*----------------------------------------------------------------------------*/
