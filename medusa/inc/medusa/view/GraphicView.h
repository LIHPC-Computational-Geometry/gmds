/*----------------------------------------------------------------------------*/
#ifndef GMDS_GRAPHICVIEW_H
#define GMDS_GRAPHICVIEW_H
/*----------------------------------------------------------------------------*/
#include <medusa/view/View.h>
/*----------------------------------------------------------------------------*/
#include <vtkSmartPointer.h>
/*----------------------------------------------------------------------------*/
#include <vector>
#include <string>
#include <vtkDataSetMapper.h>
#include <map>
#include <medusa/model/MedusaGrid.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkSelection.h>

/*----------------------------------------------------------------------------*/
class vtkRenderer;
class vtkRenderWindow;
class vtkActor;
class vtkRenderWindowInteractor;
/*----------------------------------------------------------------------------*/
namespace medusa{
    class PickingListener;
    class MouseInteractorStyle;
    class GraphicView: public View{
    public:
        enum ViewType{
            View3D, View2DX, View2DY, View2DZ
        };

        /**
         * @brief constructor. When a Graphic view is created, the picking object
         *  is created too.
         */
        GraphicView(const std::string& AName, GraphicView::ViewType  AType);
        virtual ~GraphicView();
        void selectVTKCell(int ADIm, vtkIdType ACellID);
        void refresh();
        void update();
        void setVisibleON();
        void setVisibleOFF();
        void start();

        void textMode(int AMode);

        void showSetOfTets(std::vector<vtkIdType> AIDs);

        //Visualization operations
        void wireframe();
        void solid();

        void opacity();
        void setOpacityON();
        void setOpacityOFF();


        //FF operations
        void toggleSingGraph();
        void buildSingGraph(std::vector<vtkIdType> AIDs);
        void buildSingLine(int AID, std::vector<vtkIdType> AIDs);
        void resetSingLine();

        void unlockLines();
        int selectLine(const vtkSmartPointer<vtkActor> ADeletedActor);


        //Dual operations
        void createSurface(std::vector<vtkIdType> AIDs, int ASheetID);
        void createPlanSurface(int ASheetID, MedusaGrid* ASurf);

        void viewSurface();

        void updateDual(int ANbZones);
        void setDualVisibilityON();
        void setDualVisibilityOFF();

        void unlockSheet();

        int removeActor(vtkSmartPointer<vtkActor> ADeletedActor);
        void removeActor(int ADeletedID);

        void showAxis(double ACoords[3][3]);
        void removeAxis();

        void setModeDomain();
        void setModeDual();

        void updateMesh();

        void viewBlocks();
        void toggleBlocks();


        //Geom operation
        void surfacePicked(std::vector<vtkIdType> AIDs);
        void removeLastSelection();
        void removeSelection();


        void updateCut();
        void toggleCutOpacity();

        void toggleBRep();


        void hideSetOfTets();
        void addTetSetActor();

        void showPoints(std::vector<std::vector<double>> coords);


        void resetDual();

    private:
        std::string m_name;
        bool m_visible;
        ViewType m_type;
        vtkSmartPointer<vtkRenderer> m_vtk_renderer;
        vtkSmartPointer<vtkRenderWindow> m_vtk_window;
        vtkSmartPointer<vtkRenderWindowInteractor> m_vtk_interactor;
        vtkSmartPointer<MouseInteractorStyle> m_mouse;
        std::vector<vtkSmartPointer<vtkActor> > m_vtk_actors;
        std::map<vtkSmartPointer<vtkActor>,int> m_vtk_surface_actors;
        std::map<vtkSmartPointer<vtkActor>,int> m_vtk_tets_actors;
        std::map<vtkSmartPointer<vtkActor>,int> m_vtk_sing_actors;

        vtkSmartPointer<vtkActor> m_vtk_tets_set;
        std::vector<vtkSmartPointer<vtkActor> > m_vtk_sets_actors;

        vtkSmartPointer<vtkActor> m_vtk_cut_actor;
        std::vector<vtkSmartPointer<vtkActor>> m_vtk_selection_actors;

        vtkSmartPointer<vtkDataSetMapper> selectedMapper;
        vtkSmartPointer<vtkActor> selectedActor;
        vtkDataSetMapper *m_mapper_dual;
        vtkSmartPointer<vtkActor> m_singGraph;

        vtkSmartPointer<vtkTextActor> m_textActor;

        std::map<int, std::vector<double>> m_surf_to_colors;


        void wireframeAll();

        bool m_opacity = true;

        int m_op_cut = 1;

        int m_nb_zones = 0;

        int m_nb_surfaces = 0;


        //The mode of the view:
        /* 0 = domain editing : show the domain tet mesh without data colormap + the dual surfaces with colors
         * 1 = dual structure editing : show the tet mesh with dual data color, dual zones + dual surfaces in red
         */
        int m_mode = 0;

        int surface_mode = 0;

        bool m_block_mode = false;

        /*
         * -1 not initialized, 0 hidden, 1 visible
         */
        int m_sinGraph_visibility = -1;

        int test = 0;

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_GRAPHICVIEW_H
/*----------------------------------------------------------------------------*/

