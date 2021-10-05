/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDUSAGRID_H
#define GMDS_MEDUSAGRID_H
/*----------------------------------------------------------------------------*/
#include <vtkUnstructuredGrid.h>
#include <gmds/ig/Mesh.h>
#include <map>
#include <vtkIntArray.h>
/*----------------------------------------------------------------------------*/
namespace medusa{
    class MedusaGrid{
    public:
        MedusaGrid(const std::string& AFileName);
        MedusaGrid();
        MedusaGrid(const gmds::Mesh* AMesh, int AMode);

        virtual ~MedusaGrid();

        vtkUnstructuredGrid* grid();

        gmds::TCellID getGMDSCellID(const int ADim, const vtkIdType& AVTKId);
        vtkIdType getVTKCellID(const int ADim, gmds::TCellID AGMDSId);

        gmds::Mesh* getMesh();

        std::vector<gmds::TCellID> getSingularGraph();

        void updateDualData();
        void updateCutData();

        void updateGrid(int test);

        void UpdateBRep();

        int getNbGeomSurfs();

        std::map<int, std::vector<gmds::TCellID>> getSingLines();

        std::vector<gmds::TCellID> getOffset(int ADepth, std::vector<gmds::TCellID> AIDs);

        bool checkOffsetBoundary(std::vector<gmds::TCellID> AIDs);

        std::vector<gmds::TCellID> getNodes(std::vector<gmds::TCellID> ATetIDs);

    public:
        std::string name;

    private:
        gmds::Mesh m_gmds_mesh;
        vtkUnstructuredGrid* m_vtk_grid;
        //we store the access to the mapping from vtk to gmds i-dim cell ids
        vtkIdTypeArray* m_vtk_to_gmds[4];

        vtkSmartPointer<vtkIntArray> m_vtk_blocks;
        vtkIntArray* m_vtk_cut;
        vtkIntArray* m_vtk_BND;
        std::map<gmds::TCellID,vtkIdType> m_gmds_to_vtk[4];

        int tet_treated = m_gmds_mesh.newMark<gmds::Region>();

        int m_nb_geom_surfs = 0;

        int m_mode = 0;

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDUSAGRID_H
/*----------------------------------------------------------------------------*/
