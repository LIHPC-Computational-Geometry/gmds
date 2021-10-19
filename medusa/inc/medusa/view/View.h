/*----------------------------------------------------------------------------*/
#ifndef GMDS_VIEW_H
#define GMDS_VIEW_H
/*----------------------------------------------------------------------------*/
#include <medusa/view/MedusaObserver.h>
/*----------------------------------------------------------------------------*/
#include <vtkType.h>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace medusa {
    class View: public MedusaObserver {

    public:
        View(){;}
        virtual ~View(){;}
        virtual void refresh()=0;
        virtual void selectVTKCell(int ADIm, vtkIdType ACellID)=0;
        virtual void createSurface(const std::vector<vtkIdType> AIDs, int ASheetID)=0;
    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_VIEW_H
/*----------------------------------------------------------------------------*/

