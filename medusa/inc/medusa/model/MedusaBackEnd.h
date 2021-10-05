/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDUSABACKEND_H
#define GMDS_MEDUSABACKEND_H
/*----------------------------------------------------------------------------*/
#include <vector>
#include <string>
/*----------------------------------------------------------------------------*/
#include <medusa/view/MedusaObserver.h>
#include <medusa/model/MedusaGrid.h>
#include <gmds/dualBlocking/DualBlockingSession.h>
/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/graph/TopologyGraph.h>
#include <gmds/frame3d/Params.h>
#include <gmds/cad/GeomSmoother.h>
/*----------------------------------------------------------------------------*/

namespace medusa{
    class MedusaBackend{
    public:
        /**
         * @brief initialize the unique instance of this class. It is mandatory
         *        to call it since the getInstance() method does not test the
         *        existence of this class instance.
         */
        static void initialize();

        /**
         * @brief returns the unique instance of MedusaBackend
         */
        static MedusaBackend* getInstance();

        virtual ~MedusaBackend();
        void loadVTK(const std::string& AFileName);
        void loadVTK(const std::string &AFileName, bool FF, double AAngle, std::vector<int> ASurfs);
        void addObserver(MedusaObserver* AObs);
        std::vector<MedusaGrid*> grids();

        std::vector<vtkIdType> gmdsToVtkIDs(std::vector<gmds::TCellID> AIDs);
        std::vector<gmds::TCellID> vtkToGMDSIDs(std::vector<vtkIdType> AIDs);
        std::vector<gmds::TCellID> getOffset(int ADepth, std::vector<gmds::TCellID> AIDs);

        void modeBlocks();

        //Frame operation
        std::vector<vtkIdType> getSingGraph();

        int getNbSingLines();

        bool getFrameAxis(double ACoords[][3],vtkIdType ACellID);
        void frameFieldInit(bool FF);
        void generateFrameField();

        std::vector<gmds::TCellID> pickSingLine(int ALineID);
        int unconstrainFrameTets(std::vector<gmds::TCellID> AIDs);

        //Dual operations
        std::vector<vtkIdType> createSurface(int ADim, const vtkIdType &AVTKId, int &ASheetID, int AAxis);
        std::vector<vtkIdType> createBoundary(int ADim, const vtkIdType &AVTKId, int &sheet_ID);
        void removeSheet(int ASheetID);

        int generateDual();
        void resetDual();

        void generateBlocks();
        void resetBlock();

        int reGenerateDual();
        int correctDualRegions();


        //Geom operations
        int pickSurface(gmds::TCellID AID, std::vector<vtkIdType>& AVtkIDs);
        void removeSurfFromSelection();

        int generateCut();

        void resetCut();

        void undoCut();

        void noCut();

        std::vector<std::vector<double>> wrongGeomPoint();



        std::vector<gmds::TCellID> getTetsOfSurface(int AID);

    private:
        MedusaBackend();
        void notify();
        void initGeometry(const double AAngle=0.72);
        void startBlockingSession();
        void startCuttingSession();
        void BRep();

        void buildSingLines();

        bool checkOffsetBoundary(std::vector<gmds::TCellID> AIDs);

        void frameFieldGenerateInit(gmds::ParamsGlobal *AParamGlobal,gmds::ParamsFrameField *AParamFF, gmds::ParamsMark *AParamMark);
        void freeFrameFieldMarks(gmds::ParamsMark *AParamsMark);


    private:
        static MedusaBackend* m_instance;
        std::vector<MedusaObserver*> m_observer;
        std::vector<MedusaGrid*> m_grids;
        db::DualBlockingSession* m_session = nullptr;

        gmds::graph::TopologyGraph* m_geom_cut = nullptr;

        gmds::cad::FACManager manager;
        gmds::cad::GeomMeshLinker linkerG_T;

        std::map<int, std::vector<gmds::TCellID>> m_sing_lines;
        gmds::Variable<int>* m_free_bnd;


        std::vector<gmds::TCellID> getWrongGeomPoint();

        std::vector<vtkIdType> getVtkNode(std::vector<gmds::TCellID> AIDs);

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDUSABACKEND_H
/*----------------------------------------------------------------------------*/
