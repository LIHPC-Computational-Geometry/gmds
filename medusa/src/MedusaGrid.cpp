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
    vtkReader.setDataOptions(gmds::N|gmds::R,false);
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
        m_gmds_to_vtk[3].emplace(gmds_cell_id,vtk_cell_id);


        m_vtk_grid->InsertNextCell(vtk_type,idList);

        vtk_cell_id++;
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[3]);

    m_vtk_blocks = vtkSmartPointer<vtkIntArray>::New();
    m_vtk_blocks->SetName("blocks");
    m_vtk_blocks->SetNumberOfValues(m_gmds_mesh.getNbRegions());

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
        m_vtk_blocks->InsertValue(i,0);
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_blocks);

    m_vtk_cut = vtkIntArray::New();
    m_vtk_cut->SetName("cut");
    m_vtk_cut->SetNumberOfValues(m_gmds_mesh.getNbRegions());

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
        m_vtk_cut->InsertValue(i,1);
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_cut);

}
/*----------------------------------------------------------------------------*/
void MedusaGrid::updateCutData(){
    gmds::Variable<int>* cut;
    if(m_mode == 0) {
        cut = m_gmds_mesh.getVariable<int, GMDS_REGION>("CUT_tet");
    }else if(m_mode == 1){
        cut = m_gmds_mesh.getVariable<int, GMDS_FACE>("CUT_face");
    }

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){

        TCellID gmds_id = getGMDSCellID(3,i);
        int value = cut->value(gmds_id);
        m_vtk_cut->SetValue(i,value);
    }
}
/*----------------------------------------------------------------------------*/
MedusaGrid::MedusaGrid(const Mesh* AMesh, int AMode):m_gmds_mesh(*AMesh) {

    /*
     * AMode =
     * 0 -> normal mode, building cells
     * 1 -> BRep mode, building only surf faces
     */

    m_mode = AMode;
    int nbPoints, nbCells;

    MeshDoctor doc(&m_gmds_mesh);
    //doc.buildFacesAndR2F();
    //doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    m_vtk_grid= vtkUnstructuredGrid::New();

    for (auto i = 0; i < 4; i++) {
        m_vtk_to_gmds[i] = vtkIdTypeArray::New();
        m_vtk_to_gmds[i]->SetName("GMDS_ID");
    }

    if (m_mode == 0) {

        nbPoints = m_gmds_mesh.getNbNodes();
        nbCells = m_gmds_mesh.getNbRegions();


        //=============== POINTS CREATION ===============
        vtkPoints *grid_points = vtkPoints::New();
        grid_points->SetNumberOfPoints(nbPoints);
        grid_points->SetDataTypeToDouble();
        m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

        //map to store the correspondence from a gmds node id to a vtk point id

        std::map<TCellID, int> map_node_ids;
        int vtk_node_id = 0;
        for (auto gmds_node_id:m_gmds_mesh.nodes()) {
            //id mapping
            map_node_ids[gmds_node_id] = vtk_node_id;
            //point definition
            math::Point loc = m_gmds_mesh.get<Node>(gmds_node_id).getPoint();
            grid_points->InsertPoint(vtk_node_id, loc.X(), loc.Y(), loc.Z());
            m_vtk_to_gmds[0]->InsertValue(vtk_node_id, gmds_node_id);
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
        m_gmds_to_vtk[3].emplace(gmds_cell_id, vtk_cell_id);


        m_vtk_grid->InsertNextCell(vtk_type, idList);

        vtk_cell_id++;
    }


    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[3]);
    }else if(m_mode == 1) {
        gmds::Variable<int>* BND_Color = m_gmds_mesh.getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");


        std::set<TCellID> nodes_set;
        std::vector<TCellID> BND_faces;

        for(auto f : m_gmds_mesh.faces()){
            if(BND_Color->value(f) != 0){
                Face face = m_gmds_mesh.get<Face>(f);
                BND_faces.push_back(f);
                for(auto n : face.getIDs<Node>()){
                    nodes_set.insert(n);
                }
            }
        }

        nbPoints = nodes_set.size();
        nbCells = BND_faces.size();

        //=============== POINTS CREATION ===============
        vtkPoints* grid_points = vtkPoints::New();
        grid_points->SetNumberOfPoints(nbPoints);
        grid_points->SetDataTypeToDouble();
        m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

        std::map<TCellID, int> map_node_ids;
        int vtk_node_id=0;
        for(auto gmds_node_id:nodes_set){
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

        //=============== FACES CREATION ===============

        m_vtk_to_gmds[2]->SetNumberOfValues(nbCells);

        int vtk_face_id = 0;
        for (auto gmds_face_id:BND_faces) {
            Face f = m_gmds_mesh.get<Face>(gmds_face_id);
            std::vector<TCellID> node_ids = f.getIDs<Node>();

            vtkIdList *idList = vtkIdList::New();
            idList->SetNumberOfIds(node_ids.size());
            int i_vtk_node = 0;
            for (auto i_node:node_ids) {
                idList->SetId(i_vtk_node++, map_node_ids[i_node]);
            }
            m_vtk_to_gmds[2]->InsertValue(vtk_face_id, gmds_face_id);
            m_gmds_to_vtk[2].emplace(gmds_face_id, vtk_face_id);


            m_vtk_grid->InsertNextCell(7,idList);
            vtk_face_id++;
        }

        m_vtk_BND = vtkIntArray::New();
        m_vtk_BND->SetName("BND_Color");
        m_vtk_BND->SetNumberOfValues(BND_faces.size());

        for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
            if(BND_Color->value(m_vtk_to_gmds[2]->GetValue(i)) > m_nb_geom_surfs){
                m_nb_geom_surfs = BND_Color->value(m_vtk_to_gmds[2]->GetValue(i));
            }
            m_vtk_BND->InsertValue(i,BND_Color->value(m_vtk_to_gmds[2]->GetValue(i)));
        }

        m_vtk_grid->GetCellData()->AddArray(m_vtk_BND);
        m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[2]);
    }

    m_vtk_cut = vtkIntArray::New();
    m_vtk_cut->SetName("cut");
    m_vtk_cut->SetNumberOfValues(nbCells);

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
        m_vtk_cut->InsertValue(i,1);
    }
}
/*----------------------------------------------------------------------------*/
MedusaGrid::~MedusaGrid() {
    m_vtk_grid->Delete();
    m_vtk_to_gmds[0]->Delete();
    m_vtk_to_gmds[1]->Delete();
    m_vtk_to_gmds[2]->Delete();
    m_vtk_to_gmds[3]->Delete();
    m_vtk_BND->Delete();
    m_vtk_cut->Delete();
    m_vtk_blocks->Delete();
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
getVTKCellID(const int ADim, gmds::TCellID AGMDSId)
{
    return m_gmds_to_vtk[ADim][AGMDSId];
}
/*----------------------------------------------------------------------------*/
std::vector<gmds::TCellID> MedusaGrid::getSingularGraph(){
    Variable<int>* sing;
    std::vector<gmds::TCellID> sing_graph;
    try {
        sing = m_gmds_mesh.getVariable<int,GMDS_REGION>("sing_tet");
        for(auto r : m_gmds_mesh.regions()){
            if(sing->value(r) > 0) {
                sing_graph.push_back(r);
            }
        }
    }catch (GMDSException& e){
        std::cout<<"ERROR: Singularity variable on mesh not defined."<<std::endl;
    }

    return sing_graph;
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

    std::cout<<"Size "<<m_vtk_blocks->GetSize()<<std::endl;
    std::cout<<"Max value "<<m_vtk_blocks->GetDataTypeValueMax()<<std::endl;
    std::cout<<"NbValues "<<m_vtk_blocks->GetNumberOfValues()<<std::endl;

}
/*----------------------------------------------------------------------------*/
void MedusaGrid::updateGrid(int test){

    m_vtk_grid->Reset();

    for(int i = 0; i<4; i++)
        m_gmds_to_vtk[i].clear();

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
        //=============== FACES CREATION ===============
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
            m_gmds_to_vtk[2].emplace(gmds_face_id, vtk_face_id);


            m_vtk_grid->InsertNextCell(7,idList);
            vtk_face_id++;
        }
        m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[2]);

        //std::cout << "nb VTK faces " << m_vtk_to_gmds[2]->GetNumberOfValues() << std::endl;
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
            m_gmds_to_vtk[3].emplace(gmds_cell_id, vtk_cell_id);


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
/*----------------------------------------------------------------------------*/
void MedusaGrid::UpdateBRep(){

    m_vtk_grid->Reset();

    for(int i = 0; i<4; i++)
        m_gmds_to_vtk[i].clear();

    gmds::Variable<int>* BND_Color = m_gmds_mesh.getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    gmds::Variable<int>* cut_f = m_gmds_mesh.getVariable<int, GMDS_FACE>("CUT_face");

    std::set<TCellID> nodes_set;
    std::vector<TCellID> BND_faces;

    for(auto f : m_gmds_mesh.faces()){
        if(BND_Color->value(f) != 0 && cut_f->value(f) == 0){
            Face face = m_gmds_mesh.get<Face>(f);
            BND_faces.push_back(f);
            for(auto n : face.getIDs<Node>()){
                nodes_set.insert(n);
            }
        }
    }

    auto nbPoints = nodes_set.size();
    auto nbFaces = BND_faces.size();

    //=============== POINTS CREATION ===============
    vtkPoints* grid_points = vtkPoints::New();
    grid_points->SetNumberOfPoints(nbPoints);
    grid_points->SetDataTypeToDouble();
    m_vtk_to_gmds[0]->SetNumberOfValues(nbPoints);

    std::map<TCellID, int> map_node_ids;
    int vtk_node_id=0;
    for(auto gmds_node_id:nodes_set){
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

    //=============== FACES CREATION ===============

    m_vtk_to_gmds[2]->SetNumberOfValues(nbFaces);

    int vtk_face_id = 0;
    for (auto gmds_face_id:BND_faces) {
        Face f = m_gmds_mesh.get<Face>(gmds_face_id);
        std::vector<TCellID> node_ids = f.getIDs<Node>();

        vtkIdList *idList = vtkIdList::New();
        idList->SetNumberOfIds(node_ids.size());
        int i_vtk_node = 0;
        for (auto i_node:node_ids) {
            idList->SetId(i_vtk_node++, map_node_ids[i_node]);
        }
        m_vtk_to_gmds[2]->InsertValue(vtk_face_id, gmds_face_id);
        m_gmds_to_vtk[2].emplace(gmds_face_id, vtk_face_id);


        m_vtk_grid->InsertNextCell(7,idList);
        vtk_face_id++;
    }

    m_vtk_BND = vtkIntArray::New();
    m_vtk_BND->SetName("BND_Color");
    m_vtk_BND->SetNumberOfValues(BND_faces.size());

    for(vtkIdType i = 0; i < m_vtk_grid->GetNumberOfCells(); i++){
        m_vtk_BND->InsertValue(i,BND_Color->value(m_vtk_to_gmds[2]->GetValue(i)));
    }

    m_vtk_grid->GetCellData()->AddArray(m_vtk_BND);
    m_vtk_grid->GetCellData()->AddArray(m_vtk_to_gmds[2]);

}
/*----------------------------------------------------------------------------*/
int MedusaGrid::getNbGeomSurfs() {
    return m_nb_geom_surfs;
}
/*----------------------------------------------------------------------------*/
std::map<int,std::vector<TCellID>> MedusaGrid::getSingLines(){

    Variable<int>* sing;
    std::map<int,std::vector<gmds::TCellID>> sing_Lines;
    int line_ID = 0;
    try {
        sing = m_gmds_mesh.getVariable<int,GMDS_REGION>("sing_tet");
        m_gmds_mesh.unmarkAll<Region>(tet_treated);
        for(auto r : m_gmds_mesh.regions()){
            if(sing->value(r) > 0 && !m_gmds_mesh.isMarked<Region>(r, tet_treated)) {
                std::vector<TCellID> sing_line;
                sing_line.push_back(r);
                m_gmds_mesh.mark<Region>(r, tet_treated);
                Region tet = m_gmds_mesh.get<Region>(r);
                std::vector<TCellID> next_tets;
                for(auto n : tet.get<Node>()){
                    for(auto n_r : n.getIDs<Region>()){
                        if(sing->value(n_r) > 0 && !m_gmds_mesh.isMarked<Region>(n_r, tet_treated)){
                            m_gmds_mesh.mark<Region>(n_r, tet_treated);
                            next_tets.push_back(n_r);
                        }
                    }
                }
                while(!next_tets.empty()){
                    TCellID current = next_tets.back();
                    next_tets.pop_back();
                    sing_line.push_back(current);
                    Region tet_ = m_gmds_mesh.get<Region>(current);
                    for(auto n : tet_.get<Node>()){
                        for(auto n_r : n.getIDs<Region>()){
                            if(sing->value(n_r) > 0 && !m_gmds_mesh.isMarked<Region>(n_r, tet_treated)){
                                m_gmds_mesh.mark<Region>(n_r, tet_treated);
                                next_tets.push_back(n_r);
                            }
                        }
                    }
                }
                sing_Lines.emplace(line_ID, sing_line);
                line_ID++;
            }
        }

        m_gmds_mesh.unmarkAll<Region>(tet_treated);
    }catch (GMDSException e){
        std::cout<<"ERROR: Singularity variable on mesh not defined."<<std::endl;
    }


    return sing_Lines;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaGrid::getOffset(int ADepth, std::vector<TCellID> AIDs) {

    std::vector<TCellID> result = AIDs;
    std::vector<TCellID> base = AIDs; //On part de ça pour faire le offset
    std::vector<TCellID> offset;
    m_gmds_mesh.unmarkAll<Region>(tet_treated);
    for (int i = 0; i < ADepth; i++) {
        for (auto t : base) {
            m_gmds_mesh.mark<Region>(t, tet_treated);
            Region r = m_gmds_mesh.get<Region>(t);
            for (auto n : r.get<Node>()) {
                for (auto n_t : n.getIDs<Region>()) {
                    if (!m_gmds_mesh.isMarked<Region>(n_t, tet_treated)) {
                        offset.push_back(n_t);
                        m_gmds_mesh.mark<Region>(n_t, tet_treated);
                    }
                }
            }
        }
        base = offset; //A la prochaine itération on partira du offset créé
        for(auto t : offset){
            result.push_back(t);
        }
        offset.clear();
    }


    return result;
}
/*----------------------------------------------------------------------------*/
bool MedusaGrid::checkOffsetBoundary(std::vector<TCellID> AIDs){
    for(auto id : AIDs){
        Region r = m_gmds_mesh.get<Region>(id);
        for(auto f : r.get<Face>()){
            if(f.get<Region>().size() != 1){
                int cpt_tet = 0;
                for(auto n : f.get<Node>()){
                    for(auto t : n.getIDs<Region>()){
                        if(m_gmds_mesh.isMarked<Region>(t, tet_treated)) {
                            cpt_tet++;
                            break;//On passe au noeud d'après
                        }
                    }
                }
                
            }
        }
    }

    m_gmds_mesh.unmarkAll<Region>(tet_treated);
    return true;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaGrid::getNodes(std::vector<TCellID> ATetIDs){
    std::set<TCellID> nodes;
    for(auto t : ATetIDs){
        Region r = m_gmds_mesh.get<Region>(t);
        for(auto n : r.getIDs<Node>()){
            nodes.insert(n);
        }
    }
    std::vector<TCellID> nodes_vec;
    nodes_vec.reserve(nodes.size());
    for(auto n : nodes){
        nodes_vec.push_back(n);
    }

    return nodes_vec;
}