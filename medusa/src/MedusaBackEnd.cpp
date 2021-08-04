/*----------------------------------------------------------------------------*/
#include <medusa/model/MedusaBackEnd.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>

#include <gmds/frame3d/FieldGenerator.h>
#include <gmds/frame3d/OpenNLFieldSolverStrategy.h>
#include <gmds/io/VTKWriter.h>

#include <utility>
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
    m_grids[0]->name = "Tets";
    int nb_isolated_node = 0;
    int nb_nodes = m_grids[0]->getMesh()->getNbNodes();
    std::cout<<nb_nodes<<std::endl;


    for(auto n : m_grids[0]->getMesh()->nodes()){
        Node node = m_grids[0]->getMesh()->get<Node>(n);
        if(node.nbRegions() == 0){
            std::cout<<"isolated node"<<std::endl;
            nb_isolated_node++;
            m_grids[0]->getMesh()->deleteNode(n);
        }
    }
    int nb_nodes2 = m_grids[0]->getMesh()->getNbNodes();
    std::cout<<nb_nodes2<<std::endl;
    //startBlockingSession();
    initGeometry();

    IGMeshIOService ioService_m(m_grids[0]->getMesh());
    VTKWriter vtkWriter_m(&ioService_m);
    vtkWriter_m.setCellOptions(gmds::N|gmds::R);
    vtkWriter_m.setDataOptions(gmds::N|gmds::R);
    vtkWriter_m.write("/home/simon/Data/Results_debug/test.vtk");

    frameFieldInit(false);
    buildSingLines();
    startCuttingSession();
    BRep();
    notify();


}
/*----------------------------------------------------------------------------*/
void MedusaBackend::loadVTK(const std::string &AFileName, bool FF, double AAngle, std::vector<int> ASurfs) {

    m_grids.push_back(new MedusaGrid(AFileName));
    m_grids[0]->name = "Tets";
    //startBlockingSession();
    initGeometry(AAngle);
    if(!ASurfs.empty()){
        for(auto s : ASurfs)
            unconstrainFrameTets(getTetsOfSurface(s));
    }
    frameFieldInit(FF);
    buildSingLines();
    startCuttingSession();
    BRep();
    notify();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::modeBlocks(){
    if(m_grids.size() != 1) {
        std::cout<<"grid name = "<<m_grids[1]->name<<std::endl;
        m_grids.pop_back();
    }
    startBlockingSession();
    notify();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::BRep(){
    MedusaGrid* BRep = new MedusaGrid(m_grids[0]->getMesh(), 1);
    BRep->name = "BRep";
    m_grids.push_back(BRep);
}

/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::createSurface(const int ADim, const vtkIdType &AVTKId, int &sheet_ID, int AAxis){

    TCellID gmds_id = m_grids[0]->getGMDSCellID(ADim,AVTKId);

    MedusaGrid* surface_grid = new MedusaGrid();
    std::vector<TCellID> surf_ids = m_session->createSurface(gmds_id,AAxis,sheet_ID,surface_grid->getMesh());
    std::vector<vtkIdType> gmds_to_vtk;
    if(sheet_ID != -1) {
        for (auto id : surf_ids) {
            gmds_to_vtk.push_back(m_grids[0]->getVTKCellID(3, id));
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
        gmds_to_vtk.push_back(m_grids[0]->getVTKCellID(3, id));
    }
    surface_grid->updateGrid(2);
    m_grids.push_back(surface_grid);
    return gmds_to_vtk;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::startBlockingSession(){
    MedusaGrid* block = new MedusaGrid();
    block->name = "blocks";
    block->getMesh()->changeModel(MeshModel(DIM3|R|F|E|N|R2N|F2N | E2N |
                                   R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));
    m_grids.push_back(block);
    m_session = new db::DualBlockingSession(m_grids[0]->getMesh(), block->getMesh(),manager,linkerG_T);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::startCuttingSession() {
    m_geom_cut = new graph::TopologyGraph(&manager,&linkerG_T,m_grids[0]->getMesh());
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::generateCut(){
    if(m_geom_cut->generateCut() == 0) {
        m_grids[0]->updateCutData();
        m_grids[1]->UpdateBRep();
        return 0;
    }
    return -1;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetCut(){
    m_geom_cut->resetCut();
    m_grids[0]->updateCutData();
    m_grids[1]->UpdateBRep();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::undoCut(){
    m_geom_cut->undoLastCut();
    m_grids[0]->updateCutData();
    m_grids[1]->UpdateBRep();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::removeSurfFromSelection(){
    m_geom_cut->undoSelection();
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::pickSurface(TCellID AID, std::vector<vtkIdType> &AVtkIDs){
    std::vector<TCellID> ids;
    int surf_id = m_geom_cut->pickSurf(AID, ids);
    if(surf_id == -1) return -1;
    for(auto id : ids){
        AVtkIDs.push_back(m_grids[0]->getVTKCellID(3, id));
    }
    return surf_id;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaBackend::getTetsOfSurface(int AID){
    std::vector<TCellID> ids;
    m_geom_cut->pickSurf(AID, ids);
    return ids;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::removeSheet(int ASheetID){
    m_session->deleteSheet(ASheetID);
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::generateDual(){
    int nbZones=0;
    m_session->refineSurfaceSheet();
    m_session->testLissageSheet();
    m_grids[0]->updateGrid(0);
    if(m_session->colorDual(nbZones)) {
        m_grids[0]->updateDualData();
    }else{
        nbZones = -2;
    }
    m_grids[0]->updateDualData();

    return nbZones;
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::reGenerateDual(){
    int nbZones = m_session->recolorDual();
    m_grids[0]->updateDualData();
    return nbZones;
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::correctDualRegions(){
    return m_session->removeWrongDualRegions();
}
/*----------------------------------------------------------------------------*/
bool MedusaBackend::getFrameAxis(double ACoords[][3],vtkIdType ACellID){
    TCellID gmds_id = m_grids[0]->getGMDSCellID(3,ACellID);

    return m_session->getFrameAxis(ACoords,gmds_id);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::generateBlocks(){
    m_session->createBlock();
    m_grids[1]->updateGrid(1);
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetDual() {
    m_session->hardResetDual();
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::resetBlock() {
    m_session->resetBlocks();
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::getSingGraph() {
    std::vector<TCellID> gmds_graph = m_grids[0]->getSingularGraph();
    std::vector<vtkIdType> vtk_graph;
    for(auto g_id : gmds_graph){
        vtk_graph.push_back(m_grids[0]->getVTKCellID(3, g_id));
    }
    return vtk_graph;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::initGeometry(double AAngle){
    manager.initAndLinkFrom3DMesh(m_grids[0]->getMesh(),&linkerG_T, AAngle);

    IGMeshIOService ioService_m(&manager.getMeshView());
    VTKWriter vtkWriter_m(&ioService_m);
    vtkWriter_m.setCellOptions(gmds::N|gmds::F);
    vtkWriter_m.setDataOptions(gmds::N|gmds::F);
    vtkWriter_m.write("/home/simon/Data/Results_debug/BRep.vtk");
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::frameFieldInit(bool FF){
    m_free_bnd = m_grids[0]->getMesh()->newVariable<int,GMDS_NODE>("free_bnd_nodes");

    if(!FF)
        generateFrameField();


    Variable<math::Chart>* vertex_chart = m_grids[0]->getMesh()->getOrCreateVariable<math::Chart, GMDS_NODE>("chart_field");

    Variable<math::Vector3d>* var_X_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_X_POS");
    Variable<math::Vector3d>* var_Y_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_Y_POS");
    Variable<math::Vector3d>* var_Z_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_Z_POS");

    for(auto n : m_grids[0]->getMesh()->nodes()){
        math::Vector3d x = var_X_field->value(n);
        math::Vector3d y = var_Y_field->value(n);
        math::Vector3d z = var_Z_field->value(n);
        math::Chart c(x,y,z);
        vertex_chart->set(n, c);
    }
    }
/*----------------------------------------------------------------------------*/
void MedusaBackend::generateFrameField() {

    //Ici on build le frame field sur le mesh en input
    ParamsGlobal pg;
    ParamsFrameField pf;
    ParamsMark pm;

    frameFieldGenerateInit(&pg,&pf,&pm);

    OpenNLFieldSolverStrategy* solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver, m_grids[0]->getMesh(), pg, pf, pm);
    ffg.relaxBoundary(m_free_bnd);
    ffg.execute();

    buildSingLines();

    freeFrameFieldMarks(&pm);

    Variable<math::Chart>* vertex_chart = m_grids[0]->getMesh()->getOrCreateVariable<math::Chart, GMDS_NODE>("chart_field");

    Variable<math::Vector3d>* var_X_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_X_POS");
    Variable<math::Vector3d>* var_Y_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_Y_POS");
    Variable<math::Vector3d>* var_Z_field = m_grids[0]->getMesh()->getVariable<math::Vector3d,GMDS_NODE>("FF_Z_POS");

    for(auto n : m_grids[0]->getMesh()->nodes()){
        math::Vector3d x = var_X_field->value(n);
        math::Vector3d y = var_Y_field->value(n);
        math::Vector3d z = var_Z_field->value(n);
        math::Chart c(x,y,z);
        vertex_chart->set(n, c);
    }
    IGMeshIOService ioService_m(m_grids[0]->getMesh());
    VTKWriter vtkWriter_m(&ioService_m);
    vtkWriter_m.setCellOptions(gmds::N|gmds::R);
    vtkWriter_m.setDataOptions(gmds::N|gmds::R);
    vtkWriter_m.write("/home/simon/Data/Results_debug/result_FF.vtk");


    std::cout<<"\n\tFrame field generated"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::frameFieldGenerateInit(ParamsGlobal* AParamGlobal,
                                   ParamsFrameField* AParamFF,
                                   ParamsMark* AParamMark){

    AParamGlobal->algo_choice=ParamsGlobal::STABLE_HEX_GENERATION;
    AParamGlobal->start_from=ParamsGlobal::FF_GEN;
    AParamGlobal->stop_at=ParamsGlobal::FF_SMOOTH;
    AParamGlobal->with_debug_files=true;
    AParamGlobal->with_quad_bnd_constraint=false;

    AParamFF->solver_type=ParamsFrameField::OPENNL;
    AParamFF->with_cotangent_weights=false;
    AParamFF->with_smoothing=false;
    AParamFF->smoothing_algo=ParamsFrameField::RAY;
    AParamFF->smoothing_nb_iter=10;
    AParamFF->smoothing_epsilon=1e-6;
    AParamFF->with_mesh_adaptation=false;
    AParamFF->premeshing_sing_lines=false;

    AParamMark->mark_node_on_surf = m_grids[0]->getMesh()->newMark<Node>();
    AParamMark->mark_node_on_curv = m_grids[0]->getMesh()->newMark<Node>();
    AParamMark->mark_node_on_pnt  = m_grids[0]->getMesh()->newMark<Node>();
    AParamMark->mark_node_isolated= m_grids[0]->getMesh()->newMark<Node>();

    AParamMark->mark_node_frame   = m_grids[0]->getMesh()->newMark<Node>();
    AParamMark->mark_node_hard    = m_grids[0]->getMesh()->newMark<Node>();

    AParamMark->mark_edge_on_surf = m_grids[0]->getMesh()->newMark<Edge>();
    AParamMark->mark_edge_on_curv = m_grids[0]->getMesh()->newMark<Edge>();

    AParamMark->mark_face_on_surf = m_grids[0]->getMesh()->newMark<Face>();

    BoundaryOperator boundaryOp(m_grids[0]->getMesh());

    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

//==================================================================
// Mark boundary cells
    boundaryOp.markCellOnGeometry(AParamMark->mark_face_on_surf,
                                  AParamMark->mark_edge_on_surf,
                                  AParamMark->mark_node_on_surf,
                                  AParamMark->mark_edge_on_curv,
                                  AParamMark->mark_node_on_curv,
                                  AParamMark->mark_node_on_pnt,
                                  AParamMark->mark_node_isolated);

//==================================================================
// COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS

// Color names used to name the mesh variables are defined in
// the boundary operator class
    boundaryOp.colorEdges(AParamMark->mark_edge_on_curv, AParamMark->mark_node_on_pnt);

    boundaryOp.colorNodes(AParamMark->mark_node_on_pnt);

    for (auto n_id : m_grids[0]->getMesh()->nodes()) {
        Node n = m_grids[0]->getMesh()->get<Node>(n_id);
        if (m_grids[0]->getMesh()->isMarked(n, AParamMark->mark_node_on_surf)) {
//     math::Vector3d nv = boundaryOp.getOutputNormalOfABoundaryNode(n);
//   m_bnd_normals[n.id()] = math::Vector3d(nv.X(), nv.Y(), nv.Z());
        }
    }

// now we color nodes on curves and surfaces
    Variable<int>* color_f = m_grids[0]->getMesh()->getVariable<int, GMDS_FACE>("BND_SURFACE_COLOR");
    Variable<int>* color_c = m_grids[0]->getMesh()->getVariable<int, GMDS_EDGE>("BND_CURVE_COLOR");

    Variable<int>* color_nf = 0;
    try {
        color_nf = m_grids[0]->getMesh()->newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");
    } catch (GMDSException& e) {
        color_nf = m_grids[0]->getMesh()->getVariable<int,GMDS_NODE>("BND_SURFACE_COLOR");
    }
    Variable<int>* color_nc = 0;
    try {
        color_nc = m_grids[0]->getMesh()->newVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
    } catch (GMDSException& e) {
        color_nc = m_grids[0]->getMesh()->getVariable<int,GMDS_NODE>("BND_CURVE_COLOR");
    }

    for (auto f_id: m_grids[0]->getMesh()->faces()) {
        Face f = m_grids[0]->getMesh()->get<Face>(f_id);

// only faces on surface are of interest
        if (!m_grids[0]->getMesh()->isMarked(f, AParamMark->mark_face_on_surf))
            continue;

        std::vector<Node> f_nodes = f.get<Node>();
        for (auto ni : f_nodes) {
            if (m_grids[0]->getMesh()->isMarked(ni, AParamMark->mark_node_on_surf) &&
                !m_grids[0]->getMesh()->isMarked(ni, AParamMark->mark_node_on_curv) &&
                !m_grids[0]->getMesh()->isMarked(ni, AParamMark->mark_node_on_pnt)) {
                (*color_nf)[ni.id()] = (*color_f)[f.id()];
            }
        }
    }
    for (auto e_id:m_grids[0]->getMesh()->edges()){
        Edge e = m_grids[0]->getMesh()->get<Edge>(e_id);
// only edges on surface are of interest
        if (!m_grids[0]->getMesh()->isMarked(e, AParamMark->mark_edge_on_curv))
            continue;

        std::vector<Node> e_nodes = e.get<Node>();
        for (auto ni : e_nodes) {
            if (m_grids[0]->getMesh()->isMarked(ni, AParamMark->mark_node_on_curv) &&
                !m_grids[0]->getMesh()->isMarked(ni, AParamMark->mark_node_on_pnt)) {
                (*color_nc)[ni.id()] = (*color_c)[e.id()];
            }
        }
    }

    Variable<double>* cot_w;
    try {
        cot_w = m_grids[0]->getMesh()->getVariable<double,GMDS_EDGE>("cot_weight");
    }catch (GMDSException &e){
        cot_w = m_grids[0]->getMesh()->newVariable<double,GMDS_EDGE>("cot_weight");
    }

// All weights equal to 1
    for (auto e_id:m_grids[0]->getMesh()->edges()){

        (*cot_w)[e_id] = 1;
    }

}
/*----------------------------------------------------------------------------*/
void MedusaBackend::freeFrameFieldMarks(ParamsMark* AParamsMark) {

    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_on_surf );
    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_on_curv );
    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_on_pnt  );
    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_isolated);
    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_frame   );
    m_grids[0]->getMesh()->unmarkAll<Node>(AParamsMark->mark_node_hard    );
    m_grids[0]->getMesh()->unmarkAll<Edge>(AParamsMark->mark_edge_on_surf );
    m_grids[0]->getMesh()->unmarkAll<Edge>(AParamsMark->mark_edge_on_curv );
    m_grids[0]->getMesh()->unmarkAll<Face>(AParamsMark->mark_face_on_surf );

    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_on_surf );
    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_on_curv );
    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_on_pnt  );
    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_isolated);
    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_frame   );
    m_grids[0]->getMesh()->freeMark<Node>(AParamsMark->mark_node_hard    );
    m_grids[0]->getMesh()->freeMark<Edge>(AParamsMark->mark_edge_on_surf );
    m_grids[0]->getMesh()->freeMark<Edge>(AParamsMark->mark_edge_on_curv );
    m_grids[0]->getMesh()->freeMark<Face>(AParamsMark->mark_face_on_surf );
}
/*----------------------------------------------------------------------------*/
void MedusaBackend::buildSingLines() {
    m_sing_lines = m_grids[0]->getSingLines();
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaBackend::pickSingLine(int ALineID){
    return m_sing_lines[ALineID];
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::getNbSingLines(){
    return m_sing_lines.size();
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::gmdsToVtkIDs(std::vector<TCellID> AIDs){

    std::vector<vtkIdType> result;
    result.reserve(AIDs.size());

    for(auto id : AIDs){
        result.push_back(m_grids[0]->getVTKCellID(3, id));
    }

    return result;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaBackend::vtkToGMDSIDs(std::vector<vtkIdType> AIDs){

    std::vector<TCellID> result;
    result.reserve(AIDs.size());

    for(auto id : AIDs){
        result.push_back(m_grids[0]->getGMDSCellID(3,id));
    }

    return result;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaBackend::getOffset(int ADepth, std::vector<TCellID> AIDs){
    std::vector<TCellID> result = grids()[0]->getOffset(ADepth, AIDs);
    return result;
}
/*----------------------------------------------------------------------------*/
int MedusaBackend::unconstrainFrameTets(const std::vector<TCellID> AIDs){

    for(auto t : AIDs){
        Region r = grids()[0]->getMesh()->get<Region>(t);
        for(auto n : r.getIDs<Node>())
        m_free_bnd->set(n,1);
    }

    return 0;
}
/*----------------------------------------------------------------------------*/
bool MedusaBackend::checkOffsetBoundary(std::vector<TCellID> AIDs){

}
/*----------------------------------------------------------------------------*/
void MedusaBackend::noCut() {
    m_geom_cut->noCut();
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MedusaBackend::getWrongGeomPoint(){

    std::vector<int> wrong_points;

    std::vector<TCellID> gmds_node;

    Variable<math::Chart>* vertex_chart = m_grids[0]->getMesh()->getVariable<math::Chart, GMDS_NODE>("chart_field");

    std::vector<cad::GeomPoint*> points;
    manager.getPoints(points);

    for(auto p : points){
        for(auto s : p->surfaces()){

            std::vector<cad::GeomCurve*> curves_on_surface;
            for(auto c : p->curves()){
                if(std::find(c->surfaces().begin(),c->surfaces().end(), s) != c->surfaces().end()){
                    curves_on_surface.push_back(c);
                }
            }

            // Les courbes geom incidentes au point sur la surface
            if(curves_on_surface.size() == 2) {
                cad::GeomCurve *c1 = curves_on_surface[0];
                cad::GeomCurve *c2 = curves_on_surface[1];

                Node n;
                Node n1;
                Node n2;

                math::Vector3d face_norm;

                for (auto mesh_node : m_grids[0]->getMesh()->nodes()) {
                    if (linkerG_T.getGeomDim<Node>(mesh_node) == 1 && linkerG_T.getGeomId<Node>(mesh_node) == p->id()) {
                        n = m_grids[0]->getMesh()->get<Node>(mesh_node);
                    }
                }

                //POUR NOUVELLE METHODE
                Edge e1;
                Edge e2;

                for (auto e_n : n.get<Edge>()) {
                    std::vector<Node> node_e = e_n.get<Node>();
                    if (node_e[0].id() == n.id()) {
                        if (linkerG_T.getGeomDim(node_e[1]) == 2) {
                            if (linkerG_T.getGeomId(node_e[1]) == c1->id()) {
                                n1 = node_e[1];
                                for (auto f : e_n.get<Face>()) {
                                    if (linkerG_T.getGeomDim(f) == 3 && linkerG_T.getGeomId(f) == s->id()) {
                                        face_norm = f.normal();
                                    }
                                }
                            } else if (linkerG_T.getGeomId(node_e[1]) == c2->id()) {
                                n2 = node_e[1];
                                for (auto f : e_n.get<Face>()) {
                                    if (linkerG_T.getGeomDim(f) == 3 && linkerG_T.getGeomId(f) == s->id()) {
                                        face_norm = f.normal();
                                    }
                                }
                            }
                        }
                    } else {
                        if (linkerG_T.getGeomDim(node_e[0]) == 2) {
                            if (linkerG_T.getGeomId(node_e[0]) == c1->id()) {
                                n1 = node_e[0];
                                for (auto f : e_n.get<Face>()) {
                                    if (linkerG_T.getGeomDim(f) == 3 && linkerG_T.getGeomId(f) == s->id()) {
                                        face_norm = f.normal();
                                    }
                                }
                            } else if (linkerG_T.getGeomId(node_e[0]) == c2->id()) {
                                n2 = node_e[0];
                                for (auto f : e_n.get<Face>()) {
                                    if (linkerG_T.getGeomDim(f) == 3 && linkerG_T.getGeomId(f) == s->id()) {
                                        face_norm = f.normal();
                                    }
                                }
                            }
                        }
                    }
                }

                math::Chart chart = vertex_chart->value(n.id());
                math::Chart chart1 = vertex_chart->value(n1.id());
                math::Chart chart2 = vertex_chart->value(n2.id());

                double max = 0;
                int i_max = 0;
                for (int i = 0; i < 3; i++) {
                    math::Vector3d vec = chart1[i];
                    double product = vec.dot(face_norm);
                    product = fabs(product);
                    if (product > max) {
                        max = product;
                        i_max = i;
                    }
                }


                //----------------------------------------------------------------------------------------------
                //TEST D'UNE NOUVELLE METHODE






                //----------------------------------------------------------------------------------------------

                std::vector<math::Vector3d> starting_vecs;
                if (i_max == 0) {
                    starting_vecs.push_back(chart1[1]);
                    starting_vecs.push_back(chart1[2]);
                } else if (i_max == 1) {
                    starting_vecs.push_back(chart1[0]);
                    starting_vecs.push_back(chart1[2]);
                } else {
                    starting_vecs.push_back(chart1[0]);
                    starting_vecs.push_back(chart1[1]);
                }

                bool wrong_point = false;
                for (auto starting_vec : starting_vecs) {

                    //On propage cette direction correspondante entre n2 et n0
                    math::Vector3d closest_vec_1_to_0;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart[i];
                        double product = vec.dot(starting_vec);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_1_to_0 = vec;
                        }
                    }
                    if (closest_vec_1_to_0.dot(starting_vec) < 0)
                        closest_vec_1_to_0 = -1 * closest_vec_1_to_0;

                    //std::cout<<"starting vec : "<<starting_vec<<std::endl;
                    //std::cout<<"closest_vec_1_0 : "<<closest_vec_1_to_0<<std::endl;

                    math::Vector3d closest_vec_0_to_2;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart2[i];
                        double product = vec.dot(closest_vec_1_to_0);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_0_to_2 = vec;
                        }
                    }
                    if (closest_vec_0_to_2.dot(closest_vec_1_to_0) < 0)
                        closest_vec_0_to_2 = -1 * closest_vec_0_to_2;


                    //std::cout<<"closest_vec_1_0 : "<<closest_vec_1_to_0<<std::endl;
                    //std::cout<<"closest_vec_0_2 : "<<closest_vec_0_to_2<<std::endl;

                    //On cherche la direction correspondante au vector 0 du frame de n1 pour n2
                    math::Vector3d closest_vec_2_to_1;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart1[i];
                        double product = vec.dot(closest_vec_0_to_2);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_2_to_1 = vec;
                        }
                    }
                    if (closest_vec_2_to_1.dot(closest_vec_0_to_2) < 0)
                        closest_vec_2_to_1 = -1 * closest_vec_2_to_1;

                    if (starting_vec != closest_vec_2_to_1) {
                        std::cout<<"starting vec : "<<starting_vec<<std::endl;
                        std::cout<<"closest_vec : "<<closest_vec_2_to_1<<std::endl;
                        //gmds_node.push_back(n.id());
                        wrong_point = true;
                        std::cout << "Wrong point "<<p->id()<<" on surf "<<s->id()<< std::endl;
                    }

                    /*
                    //On cherche la direction correspondante au vector 0 du frame de n1 pour n2
                    math::Vector3d closest_vec_1_to_2;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart2[i];
                        double product = vec.dot(starting_vec);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_1_to_2 = vec;
                        }
                    }
                    if (closest_vec_1_to_2.dot(starting_vec) < 0)
                        closest_vec_1_to_2 = -1 * closest_vec_1_to_2;

                    std::cout << closest_vec_1_to_2 << std::endl;

                    //On propage cette direction correspondante entre n2 et n0
                    math::Vector3d closest_vec_2_to_0;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart[i];
                        double product = vec.dot(closest_vec_1_to_2);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_2_to_0 = vec;
                        }
                    }
                    if (closest_vec_2_to_0.dot(closest_vec_1_to_2) < 0)
                        closest_vec_2_to_0 = -1 * closest_vec_2_to_0;

                    std::cout << closest_vec_2_to_0 << std::endl;

                    math::Vector3d closest_vec_0_to_1;
                    max = 0;
                    for (int i = 0; i < 3; i++) {
                        math::Vector3d vec = chart1[i];
                        double product = vec.dot(closest_vec_2_to_0);
                        product = fabs(product);
                        if (product > max) {
                            max = product;
                            closest_vec_0_to_1 = vec;
                        }
                    }
                    if (closest_vec_0_to_1.dot(closest_vec_2_to_0) < 0)
                        closest_vec_0_to_1 = -1 * closest_vec_0_to_1;

                    std::cout << closest_vec_0_to_1 << std::endl;

                    if (starting_vec != closest_vec_0_to_1) {
                        //gmds_node.push_back(n.id());
                        wrong_point = true;
                        std::cout << "Wrong point" << std::endl;
                    }*/


                    /*if(s->angle(p)<(1.5708/2)){
                        wrong_points.push_back(p->id());
                        break;
                    }*/
                }
                if (wrong_point) gmds_node.push_back(n.id());
            }
        }
    }


    /*for(auto n : m_grids[0]->getMesh()->nodes()){
        if(linkerG_T.getGeomDim<Node>(n) == 1) {
            if (std::find(wrong_points.begin(), wrong_points.end(), linkerG_T.getGeomId<Node>(n)) !=
                wrong_points.end()) {
                gmds_node.push_back(n);
            }
        }
    }*/

    return gmds_node;
}
/*----------------------------------------------------------------------------*/
std::vector<vtkIdType> MedusaBackend::getVtkNode(std::vector<TCellID> AIDs){
    std::vector<vtkIdType> vtk_ids;
    vtk_ids.reserve(AIDs.size());
    for(auto id : AIDs) {
        vtk_ids.push_back(m_grids[0]->getVTKCellID(0, id));
    }
    return vtk_ids;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<double>> MedusaBackend::wrongGeomPoint() {

    std::vector<std::vector<double>> result;

    std::vector<TCellID> points = getWrongGeomPoint();

    for(auto p : points){
        std::vector<double> tab;
        Node n = m_grids[0]->getMesh()->get<Node>(p);
        tab.push_back(n.getPoint().X());
        tab.push_back(n.getPoint().Y());
        tab.push_back(n.getPoint().Z());
        result.push_back(tab);
    }

    return result;
}