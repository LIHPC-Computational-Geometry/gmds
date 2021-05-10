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
        void addObserver(MedusaObserver* AObs);
        std::vector<MedusaGrid*> grids();

        std::vector<vtkIdType> createSurface(const int ADim, const vtkIdType &AVTKId, int &ASheetID, int AAxis);
        std::vector<vtkIdType> createBoundary(const int ADim, const vtkIdType &AVTKId, int &sheet_ID);
        void startBlockingSession();

        int generateDual();

        void removeSheet(int ASheetID);

        void getFrameAxis(double ACoords[][3],vtkIdType ACellID);

        void generateBlocks();

        void resetDual();

        void resetBlock();

        void smoothBlocks();
        void resetBlockSmoothing();
        std::vector<vtkIdType> getSingGraph();


    private:
        MedusaBackend();
        void notify();

    private:
        static MedusaBackend* m_instance;
        std::vector<MedusaObserver*> m_observer;
        std::vector<MedusaGrid*> m_grids;
        db::DualBlockingSession* m_session = NULL;

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDUSABACKEND_H
/*----------------------------------------------------------------------------*/
