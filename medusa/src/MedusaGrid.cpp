/*----------------------------------------------------------------------------*/
#include "medusa/model/MedusaGrid.h"
/*----------------------------------------------------------------------------*/
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
/*----------------------------------------------------------------------------*/
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <gmds/ig/MeshDoctor.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
using namespace gmds;

MedusaGrid::MedusaGrid()
        :m_gmds_mesh(MeshModel(DIM3|R|F|E|N|R2N|R2F|F2N|E2N|R2E|F2E)){

    m_vtk_grid= vtkUnstructuredGrid::New();
    for(auto i=0;i<4;i++){
        m_vtk_to_gmds[i]= vtkIdTypeArray::New();
        m_vtk_to_gmds[i]->SetName("GMDS_ID");
    }
    /*vtkPoints *pts = vtkPoints::New();

    pts->SetNumberOfPoints(8);
    pts->InsertPoint(0, 0, 0, 0);
    pts->InsertPoint(1, 1, 0, 0);
    pts->InsertPoint(2, 1, 1, 0);
    pts->InsertPoint(3, 0, 1, 0);
    pts->InsertPoint(4, 0, 0, 1);
    pts->InsertPoint(5, 1, 0, 1);
    pts->InsertPoint(6, 1, 1, 1);
    pts->InsertPoint(7, 0, 1, 1);
    vtkHexahedron* aHexahedron = vtkHexahedron::New();
    aHexahedron->GetPointIds()->SetId(0, 0);
    aHexahedron->GetPointIds()->SetId(1, 1);
    aHexahedron->GetPointIds()->SetId(2, 2);
    aHexahedron->GetPointIds()->SetId(3, 3);
    aHexahedron->GetPointIds()->SetId(4, 4);
    aHexahedron->GetPointIds()->SetId(5, 5);
    aHexahedron->GetPointIds()->SetId(6, 6);
    aHexahedron->GetPointIds()->SetId(7, 7);
    m_vtk_grid->Allocate(1, 1);
    m_vtk_grid->InsertNextCell(aHexahedron->GetCellType(),
                               aHexahedron->GetPointIds());
    m_vtk_grid->SetPoints(pts);*/
}
/*----------------------------------------------------------------------------*/
MedusaGrid::MedusaGrid(const std::string& AFileName)
        :m_gmds_mesh(MeshModel(DIM3|R|F|E|N|R2N|F2N | E2N |
                               R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E)){

    m_vtk_grid= vtkUnstructuredGrid::New();
    //================================================
    IGMeshIOService ioService(&m_gmds_mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(AFileName);

    MeshDoctor doc(&m_gmds_mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    auto nbPoints = m_gmds_mesh.getNbNodes();
    auto nbCells = m_gmds_mesh.getNbRegions();

    for(auto i=0;i<4;i++){
        m_vtk_to_gmds[i]= vtkIdTypeArray::New();
        m_vtk_to_gmds[i]->SetName("GMDS_ID");
    }



    //=============== POINTS CREATION ===============
    vtkPoints* grid_points = vtkPoints::New();
    grid_points->SetNumberOfPoints(nbPoints);
    grid_points->SetDataTypeToDouble();
    m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

    //map to store the correspondence from a gmds node id to a vtk point id

    std::map<TCellID, int> map_node_ids;
    int vtk_node_id=0;
    for(auto gmds_node_id:m_gmds_mesh.nodes()){
        //id mapping
        map_node_ids[gmds_node_id]=vtk_node_id;
        //point definition
        math::Point loc = m_gmds_mesh.get<Node>(gmds_node_id).getPoint();
        grid_points->InsertPoint(vtk_node_id,loc.X(),loc.Y(),loc.Z());
        m_vtk_to_gmds[0]->InsertValue(vtk_node_id,gmds_node_id);
        vtk_node_id++;
    }

    m_vtk_grid->SetPoints(grid_points);
    grid_points->Delete();

    m_vtk_grid->GetPointData()->AddArray(m_vtk_to_gmds[0]);


    //=============== CELLS CREATION ===============
    m_vtk_grid->Allocate(nbCells);

    m_vtk_to_gmds[3]->SetNumberOfValues(nbCells);

    int vtk_cell_id=0;
    for(auto gmds_cell_id:m_gmds_mesh.regions())
    {
        Region r = m_gmds_mesh.get<Region>(gmds_cell_id);
        int vtk_type=0;
        switch (r.type()){
            case GMDS_HEX:
                vtk_type = VTK_HEXAHEDRON;
                break;
            case GMDS_TETRA:
                vtk_type = VTK_TETRA;
                break;
            case GMDS_PRISM3:
                vtk_type = VTK_WEDGE;
                break;
            case GMDS_PYRAMID:
                vtk_type = VTK_PYRAMID;
                break;
            default:
                throw GMDSException("Input GMDS mesh contains cell not handled by VTK");
        }

        std::vector<TCellID> node_ids = r.getIDs<Node>();

        vtkIdList* idList = vtkIdList::New();
        idList->SetNumberOfIds(node_ids.size());
        int i_vtk_node = 0;
        for(auto i_node:node_ids){
            idList->SetId(i_vtk_node++, map_node_ids[i_node]);
        }
        m_vtk_to_gmds[3]->InsertValue(vtk_cell_id,gmds_cell_id);
        m_gmds_to_vtk.emplace(gmds_cell_id,vtk_cell_id);


        m_vtk_grid->InsertNextCell(vtk_type,idList);

        vtk_cell_id++;
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[3]);

    m_vtk_blocks = vtkIntArray::New();
    m_vtk_blocks->SetName("blocks");
    m_vtk_blocks->SetNumberOfValues(m_gmds_mesh.getNbRegions());

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
        m_vtk_blocks->InsertValue(i,0);
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_blocks);

}
/*----------------------------------------------------------------------------*/
MedusaGrid::MedusaGrid(const Mesh AMesh):m_gmds_mesh(AMesh){
    MeshDoctor doc(&m_gmds_mesh);
    //doc.buildFacesAndR2F();
    //doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    auto nbPoints = m_gmds_mesh.getNbNodes();
    auto nbCells = m_gmds_mesh.getNbRegions();

    for(auto i=0;i<4;i++){
        m_vtk_to_gmds[i] = vtkIdTypeArray::New();
        m_vtk_to_gmds[i]->SetName("GMDS_ID");
    }



    //=============== POINTS CREATION ===============
    vtkPoints* grid_points = vtkPoints::New();
    grid_points->SetNumberOfPoints(nbPoints);
    grid_points->SetDataTypeToDouble();
    m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

    //map to store the correspondence from a gmds node id to a vtk point id

    std::map<TCellID, int> map_node_ids;
    int vtk_node_id=0;
    for(auto gmds_node_id:m_gmds_mesh.nodes()){
        //id mapping
        map_node_ids[gmds_node_id]=vtk_node_id;
        //point definition
        math::Point loc = m_gmds_mesh.get<Node>(gmds_node_id).getPoint();
        grid_points->InsertPoint(vtk_node_id,loc.X(),loc.Y(),loc.Z());
        m_vtk_to_gmds[0]->InsertValue(vtk_node_id,gmds_node_id);
        vtk_node_id++;
    }

    m_vtk_grid->SetPoints(grid_points);
    grid_points->Delete();

    m_vtk_grid->GetPointData()->AddArray(m_vtk_to_gmds[0]);

        //=============== CELLS CREATION ===============
        m_vtk_grid->Allocate(nbCells);

        m_vtk_to_gmds[3]->SetNumberOfValues(nbCells);

        int vtk_cell_id = 0;
        for (auto gmds_cell_id:m_gmds_mesh.regions()) {
            Region r = m_gmds_mesh.get<Region>(gmds_cell_id);
            int vtk_type = 0;
            switch (r.type()) {
                case GMDS_HEX:
                    vtk_type = VTK_HEXAHEDRON;
                    break;
                case GMDS_TETRA:
                    vtk_type = VTK_TETRA;
                    break;
                case GMDS_PRISM3:
                    vtk_type = VTK_WEDGE;
                    break;
                case GMDS_PYRAMID:
                    vtk_type = VTK_PYRAMID;
                    break;
                default:
                    throw GMDSException("Input GMDS mesh contains cell not handled by VTK");
            }

            std::vector<TCellID> node_ids = r.getIDs<Node>();

            vtkIdList *idList = vtkIdList::New();
            idList->SetNumberOfIds(node_ids.size());
            int i_vtk_node = 0;
            for (auto i_node:node_ids) {
                idList->SetId(i_vtk_node++, map_node_ids[i_node]);
            }
            m_vtk_to_gmds[3]->InsertValue(vtk_cell_id, gmds_cell_id);
            m_gmds_to_vtk.emplace(gmds_cell_id, vtk_cell_id);


            m_vtk_grid->InsertNextCell(vtk_type, idList);

            vtk_cell_id++;
        }

        m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[3]);
}
/*----------------------------------------------------------------------------*/
MedusaGrid::~MedusaGrid() {
    m_vtk_grid->Delete();
    m_vtk_to_gmds[0]->Delete();
    m_vtk_to_gmds[1]->Delete();
    m_vtk_to_gmds[2]->Delete();
    m_vtk_to_gmds[3]->Delete();
}
/*----------------------------------------------------------------------------*/
vtkUnstructuredGrid* MedusaGrid::grid() {return m_vtk_grid;}
/*----------------------------------------------------------------------------*/
gmds::TCellID MedusaGrid::
getGMDSCellID(const int ADim, const vtkIdType &AVTKId)
{
    return (m_vtk_to_gmds[ADim])->GetValue(AVTKId);
}
/*----------------------------------------------------------------------------*/
gmds::Mesh* MedusaGrid::getMesh() {return &m_gmds_mesh;}
/*----------------------------------------------------------------------------*/
vtkIdType MedusaGrid::
getVTKCellID(gmds::TCellID AGMDSId)
{
    return m_gmds_to_vtk[AGMDSId];
}
/*----------------------------------------------------------------------------*/
void MedusaGrid::
updateDualData()
{

    gmds::Variable<int>* block = m_gmds_mesh.getVariable<int, GMDS_REGION>("blocks");


    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){

        TCellID gmds_id = getGMDSCellID(3,i);
        int value = (*block)[gmds_id];
        m_vtk_blocks->SetValue(i,value);
    }

}
/*----------------------------------------------------------------------------*/
void MedusaGrid::updateGrid(int test){

    m_vtk_grid->Reset();

    m_gmds_to_vtk.clear();

    auto nbPoints = m_gmds_mesh.getNbNodes();
    auto nbCells = m_gmds_mesh.getNbRegions();

    //=============== POINTS CREATION ===============
    vtkPoints* grid_points = vtkPoints::New();
    grid_points->SetNumberOfPoints(nbPoints);
    grid_points->SetDataTypeToDouble();
    m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

    std::map<TCellID, int> map_node_ids;
    int vtk_node_id=0;
    for(auto gmds_node_id:m_gmds_mesh.nodes()){
        //id mapping
        map_node_ids[gmds_node_id]=vtk_node_id;
        //point definition
        math::Point loc = m_gmds_mesh.get<Node>(gmds_node_id).getPoint();
        grid_points->InsertPoint(vtk_node_id,loc.X(),loc.Y(),loc.Z());
        m_vtk_to_gmds[0]->InsertValue(vtk_node_id,gmds_node_id);
        vtk_node_id++;
    }
    m_vtk_grid->SetPoints(grid_points);
    grid_points->Delete();

    m_vtk_grid->GetPointData()->AddArray(m_vtk_to_gmds[0]);

if(test == 2){
    //=============== CELLS CREATION ===============
    auto nbFaces = m_gmds_mesh.getNbFaces();
    m_vtk_to_gmds[2]->SetNumberOfValues(nbFaces);

    int vtk_face_id = 0;
    for (auto gmds_face_id:m_gmds_mesh.faces()) {
        Face f = m_gmds_mesh.get<Face>(gmds_face_id);
        std::vector<TCellID> node_ids = f.getIDs<Node>();

        vtkIdList *idList = vtkIdList::New();
        idList->SetNumberOfIds(node_ids.size());
        int i_vtk_node = 0;
        for (auto i_node:node_ids) {
            idList->SetId(i_vtk_node++, map_node_ids[i_node]);
        }
        m_vtk_to_gmds[2]->InsertValue(vtk_face_id, gmds_face_id);
        m_gmds_to_vtk.emplace(gmds_face_id, vtk_face_id);


        m_vtk_grid->InsertNextCell(7,idList);
        vtk_face_id++;
    }
    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[2]);

    std::cout << "nb VTK faces " << m_vtk_to_gmds[2]->GetNumberOfValues() << std::endl;
} else {
    //=============== CELLS CREATION ===============
    m_vtk_grid->Allocate(nbCells);

    m_vtk_to_gmds[3]->SetNumberOfValues(nbCells);

    int vtk_cell_id = 0;
    for (auto gmds_cell_id:m_gmds_mesh.regions()) {
        Region r = m_gmds_mesh.get<Region>(gmds_cell_id);
        int vtk_type = 0;
        switch (r.type()) {
            case GMDS_HEX:
                vtk_type = VTK_HEXAHEDRON;
                break;
            case GMDS_TETRA:
                vtk_type = VTK_TETRA;
                break;
            case GMDS_PRISM3:
                vtk_type = VTK_WEDGE;
                break;
            case GMDS_PYRAMID:
                vtk_type = VTK_PYRAMID;
                break;
            default:
                throw GMDSException("Input GMDS mesh contains cell not handled by VTK");
        }

        std::vector<TCellID> node_ids = r.getIDs<Node>();

        vtkIdList *idList = vtkIdList::New();
        idList->SetNumberOfIds(node_ids.size());
        int i_vtk_node = 0;
        for (auto i_node:node_ids) {
            idList->SetId(i_vtk_node++, map_node_ids[i_node]);
        }
        m_vtk_to_gmds[3]->InsertValue(vtk_cell_id, gmds_cell_id);
        m_gmds_to_vtk.emplace(gmds_cell_id, vtk_cell_id);


        m_vtk_grid->InsertNextCell(vtk_type, idList);

        vtk_cell_id++;
    }

    std::cout << "nb VTK cells " << m_vtk_to_gmds[3]->GetNumberOfValues() << std::endl;

    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[3]);
}
if(test == 0) {
    m_vtk_blocks->SetNumberOfValues(m_gmds_mesh.getNbRegions());
}
}
