/*----------------------------------------------------------------------------*/
#include <medusa/model/MedusaBackEnd.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>

/*----------------------------------------------------------------------------*/
using namespace medusa;
using namespace gmds;
MedusaBackend* MedusaBackend::m_instance = NULL;
/*----------------------------------------------------------------------------*/
MedusaBackend* MedusaBackend::getInstance() {return m_instance;}
/*----------------------------------------------------------------------------*/
void MedusaBackend::initialize()
{
    if(m_instance==NULL)
        m_instance = new MedusaBackend();
}
/*----------------------------------------------------------------------------*/
MedusaBackend::MedusaBackend()
{;}
/*----------------------------------------------------------------------------*/
MedusaBackend::~MedusaBackend()
{
    for(auto g: m_grids){
        delete g;
    }

}
/*----------------------------------------------------------------------------*/
std::vector<MedusaGrid*> MedusaBackend::grids() {
    return m_grids;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::addObserver(medusa::MedusaObserver *AObs) {
    m_observer.push_back(AObs);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::notify() {
    for(auto obs:m_observer){
        obs->update();
    }
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::loadVTK(const std::string &AFileName) {

    m_grids.push_back(new MedusaGrid(AFileName));
    startBlockingSession();
    notify();
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::createSurface(const int ADim, const vtkIdType &AVTKId, int &sheet_ID, int AAxis){

    TCellID gmds_id = m_grids[0]->getGMDSCellID(ADim,AVTKId);

    MedusaGrid* surface_grid = new MedusaGrid();
    std::vector<TCellID> surf_ids = m_session->createSurface(gmds_id,AAxis,sheet_ID,surface_grid->getMesh());
    std::vector<vtkIdType> gmds_to_vtk;
    if(sheet_ID != -1) {
        for (auto id : surf_ids) {
            gmds_to_vtk.push_back(m_grids[0]->getVTKCellID(id));
        }
        surface_grid->updateGrid(2);
        m_grids.push_back(surface_grid);
    }
    return gmds_to_vtk;
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::createBoundary(const int ADim, const vtkIdType &AVTKId, int &sheet_ID){

    TCellID gmds_id = m_grids[0]->getGMDSCellID(ADim,AVTKId);

    MedusaGrid* surface_grid = new MedusaGrid();
    std::vector<TCellID> surf_ids = m_session->createBoundary(gmds_id,sheet_ID,surface_grid->getMesh());
    std::vector<vtkIdType> gmds_to_vtk;
    for(auto id : surf_ids){
        gmds_to_vtk.push_back(m_grids[0]->getVTKCellID(id));
    }
    surface_grid->updateGrid(2);
    m_grids.push_back(surface_grid);
    return gmds_to_vtk;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::startBlockingSession(){
    MedusaGrid* block = new MedusaGrid();
    m_grids.push_back(block);
    m_session = new db::DualBlockingSession(m_grids[0]->getMesh(), block->getMesh());
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::removeSheet(int ASheetID){
    m_session->deleteSheet(ASheetID);
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::generateDual(){
    int nbZones=0;
    m_session->refineSurfaceSheet();
    m_grids[0]->updateGrid(0);
    m_grids[0]->updateDualData();
    m_session->colorDual(nbZones);
    m_grids[0]->updateDualData();
    return nbZones;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::getFrameAxis(double ACoords[][3],vtkIdType ACellID){
    TCellID gmds_id = m_grids[0]->getGMDSCellID(3,ACellID);
    m_session->getFrameAxis(ACoords,gmds_id);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::generateBlocks(){
    m_session->createBlock();
    m_grids[1]->updateGrid(1);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetDual() {
    m_session->resetDual();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetBlock() {
    m_session->resetBlocks();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::smoothBlocks() {
    m_session->smoothBlocks();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetBlockSmoothing() {
    m_session->resetBlockSmoothing();
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::getSingGraph() {
    std::vector<TCellID> gmds_graph = m_session->getSingularGraph();
    std::vector<vtkIdType> vtk_graph;
    for(auto g_id : gmds_graph){
        vtk_graph.push_back(m_grids[0]->getVTKCellID(g_id));
    }
    return vtk_graph;
}
