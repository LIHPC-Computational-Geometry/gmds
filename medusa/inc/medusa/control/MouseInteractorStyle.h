/*----------------------------------------------------------------------------*/
#ifndef GMDS_MOUSEINTERACTORSTYLE_H
#define GMDS_MOUSEINTERACTORSTYLE_H
/*----------------------------------------------------------------------------*/
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellPicker.h>
#include <vtkCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkNamedColors.h>
#include <vtkObjectFactory.h>
#include <vtkPlaneSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
/*----------------------------------------------------------------------------*/
namespace medusa {
    class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera {
    public:
        static MouseInteractorStyle *New();

        virtual void OnChar() override ;
        MouseInteractorStyle();
        void pick();
        void setPolicy();

        void wireframe();

        void solid();

        void opacity();

        void help();

        void generateDual();

        void viewDual();

        void remove();

        void generateBlocks();

        void surfaceMode();

        void finalMesh();

        void undo();

        void singularityGraph();
    };

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MOUSEINTERACTORSTYLE_H
/*----------------------------------------------------------------------------*/
