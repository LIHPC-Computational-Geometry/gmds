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
        MedusaGrid(gmds::Mesh AMesh);

        virtual ~MedusaGrid();

        vtkUnstructuredGrid* grid();

        gmds::TCellID getGMDSCellID(const int ADim, const vtkIdType& AVTKId);
        vtkIdType getVTKCellID(gmds::TCellID AGMDSId);

        gmds::Mesh* getMesh();
        void updateDualData();

        void updateGrid(int test);

    private:
        gmds::Mesh m_gmds_mesh;
        vtkUnstructuredGrid* m_vtk_grid;
        //we store the access to the mapping from vtk to gmds i-dim cell ids
        vtkIdTypeArray* m_vtk_to_gmds[4];

        vtkIntArray* m_vtk_blocks;
        std::map<gmds::TCellID,vtkIdType> m_gmds_to_vtk;

    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDUSAGRID_H
/*----------------------------------------------------------------------------*/
