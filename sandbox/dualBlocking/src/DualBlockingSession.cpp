//
// Created by calderans on 06/03/19.
//

#include <Predicates_psm.h>
#include <ctime>
#include <gmds/cad/GeomSurface.h>
#include <gmds/dualBlocking/DualBlockingSession.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/math/Line.h>
#include <gmds/math/Orientation.h>
#include <gmds/smoothy/LaplacianSmoother3C.h>

using namespace db;
using namespace gmds;

int sheet_ID = 1;

DualBlockingSession::DualBlockingSession(Mesh *AMesh, Mesh *AHMesh, cad::FACManager &AManager, cad::GeomMeshLinker &ALinker):m_mesh(AMesh), m_hmesh(AHMesh),manager(AManager), linkerG_T(ALinker) {


    GEO::PCK::initialize();


    Variable<math::Vector3d> *var_X_field;
    try {
        var_X_field = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_rep_X");
    }catch (GMDSException e){
        var_X_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_X_POS");
    }
    Variable<math::Vector3d> *var_Y_field;
    try {
        var_Y_field = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_rep_Y");
    }catch (GMDSException e){
        var_Y_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_Y_POS");
    }
    Variable<math::Vector3d> *var_Z_field;
    try {
        var_Z_field = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("quaternion_rep_Z");
    }catch (GMDSException e){
        var_Z_field = m_mesh->getVariable<math::Vector3d,GMDS_NODE>("FF_Z_POS");
    }

    try {
        m_vertex_chart = m_mesh->newVariable<math::Chart, GMDS_NODE>("chart_field");
    }catch (GMDSException e){
        m_vertex_chart = m_mesh->getVariable<math::Chart,GMDS_NODE>("chart_field");
    }

    try {
        m_tetra_chart = m_mesh->newVariable<math::Chart, GMDS_REGION>("center_chart");
    }catch (GMDSException e){
        m_tetra_chart = m_mesh->getVariable<math::Chart, GMDS_REGION>("center_chart");
    }

    try {
        m_sing = m_mesh->newVariable<int, GMDS_REGION>("sing_tet");
    }catch (GMDSException e) {
        m_sing = m_mesh->getVariable<int, GMDS_REGION>("sing_tet");
    }

    //(*m_sing)[75814] = 0;

    try {
        m_propagation_round = m_mesh->newVariable<int, GMDS_REGION>("propagation_round");
    }catch (GMDSException e) {
        m_propagation_round = m_mesh->getVariable<int, GMDS_REGION>("propagation_round");
    }

    try {
        var_info = m_mesh->newVariable<std::vector<DualSheetCreator::intersectInfo>,GMDS_EDGE>("var_info");
    }catch (GMDSException e) {
        var_info = m_mesh->getVariable<std::vector<DualSheetCreator::intersectInfo>,GMDS_EDGE>("var_info");
    }

    for(auto n : m_mesh->nodes()){
        math::Vector3d x = var_X_field->value(n);
        math::Vector3d y = var_Y_field->value(n);
        math::Vector3d z = var_Z_field->value(n);
        math::Chart c(x,y,z);
        m_vertex_chart->set(n, c);
    }


    for(auto r : m_mesh->regions()){
        Region tet = m_mesh->get<Region>(r);
        m_tetra_chart->set(r, computeChart(tet.get<Node>(), tet.center()));
        m_propagation_round->set(r, -1);
        //(*m_sing)[r] = 0;
    }

    try {
        m_sheet_X = m_mesh->newVariable<int, GMDS_REGION>("sheets_id_X");
        std::cout<<"new X variable"<<std::endl;
    }catch (GMDSException e) {
        m_sheet_X = m_mesh->getVariable<int, GMDS_REGION>("sheets_id_X");
        std::cout<<"read X variable"<<std::endl;
    }
    try {
        m_sheet_Y = m_mesh->newVariable<int, GMDS_REGION>("sheets_id_Y");
        std::cout<<"new Y variable"<<std::endl;
    }catch (GMDSException e){
        m_sheet_Y = m_mesh->getVariable<int, GMDS_REGION>("sheets_id_Y");
        std::cout<<"read Y variable"<<std::endl;
    }
    try {
        m_sheet_Z = m_mesh->newVariable<int, GMDS_REGION>("sheets_id_Z");
        std::cout<<"new Z variable"<<std::endl;
    }catch (GMDSException e){
        m_sheet_Z = m_mesh->getVariable<int, GMDS_REGION>("sheets_id_Z");
        std::cout<<"read Z variable"<<std::endl;
    }

    try {
        m_block = m_mesh->newVariable<int, GMDS_REGION>("blocks");
    }catch (GMDSException e) {
        m_block = m_mesh->getVariable<int, GMDS_REGION>("blocks");
    }



    //m_block = m_mesh->newVariable<int,GMDS_REGION>("blocks");

    m_block_id = m_hmesh->newVariable<int, GMDS_REGION>("blocks");

    //manager.initAndLinkFrom3DMesh(m_mesh,&linkerG_T);

    //linkerH_G = new cad::GeomMeshLinker();
    linkerH_G.setGeometry(&manager);
    linkerH_G.setMesh(m_hmesh);



    mark_visited = m_mesh->newMark<Region>();
    mark_ghost = m_mesh->newMark<Region>();

    wave_precedent = m_mesh->newMark<Edge>();
    face_treated = m_mesh->newMark<Face>();
    wave_tet_mark = m_mesh->newMark<Region>();


    m_X = m_hmesh->newVariable<int, GMDS_NODE>("X");
    m_Y = m_hmesh->newVariable<int, GMDS_NODE>("Y");
    m_Z = m_hmesh->newVariable<int, GMDS_NODE>("Z");

    min_length = 999999;

    for(auto e : m_mesh->edges()){
        Edge e_ = m_mesh->get<Edge>(e);
        if(min_length>e_.length()){
            min_length = e_.length();
        }
    }
    min_length = min_length*0.01;

    m_surfaceCreator = new DualSurfaceCreator(m_mesh,var_info,m_vertex_chart,m_propagation_round,
                                              m_sing,wave_precedent,face_treated,wave_tet_mark,min_length,linkerG_T,manager);


    m_boundaryCreator = new BoundarySurfaceCreator(m_mesh,var_info,m_vertex_chart,m_propagation_round,m_sing,wave_precedent,face_treated,wave_tet_mark,min_length,linkerG_T,manager);

    try {
        m_part = m_mesh->getVariable<int,GMDS_REGION>("CUT_tet");
    }catch (GMDSException e){
        m_part = m_mesh->newVariable<int,GMDS_REGION>("CUT_tet");
        for(auto r : m_mesh->regions()){
            m_part->set(r,1);
        }
    }

    TCoord min[3],max[3];
    std::vector<cad::GeomVolume*> volumes;
    manager.getVolumes(volumes);
    cad::GeomVolume* v = volumes[0];
    v->computeBoundingBox(min,max);
    distance = sqrt((min[0]-max[0])*(min[0]-max[0])+(min[1]-max[1])*(min[1]-max[1])+(min[2]-max[2])*(min[2]-max[2]));

}
/*------------------------------------------------------------------------*/

DualBlockingSession::~DualBlockingSession(){

}
/*------------------------------------------------------------------------*/

void DualBlockingSession::init() {

    std::cout<<"Dual Blocking session initialisation"<<std::endl;

    dual_sheets.clear();
    //dual_bounder.clear();

    for(auto r : m_mesh->regions()){

        // no sheet going through r in char direction X
        (*m_sheet_X)[r] = 0;
        // no sheet going through r in char direction Y
        (*m_sheet_Y)[r] = 0;
        // no sheet going through r in char direction Z
        (*m_sheet_Z)[r] = 0;

        (*m_block)[r] = 0;
    }

    sheet_ID = 0;
}
/*------------------------------------------------------------------------*/

std::vector<TCellID> DualBlockingSession::createSurface(TCellID AID, int AAxis, int &ASheetID, Mesh* ASurface_mesh) {

    ASheetID = sheet_ID;

    DualSheet dualsheet(sheet_ID,ASurface_mesh);



    //We do -1 on AAxis to fit the index range of the array
    math::Point norm((*m_tetra_chart)[AID][m_axis[AAxis-1]].X(),(*m_tetra_chart)[AID][m_axis[AAxis-1]].Y(),(*m_tetra_chart)[AID][m_axis[AAxis-1]].Z());

    m_surfaceCreator->setSheetID(sheet_ID);
    m_surfaceCreator->setTetID(AID);
    m_surfaceCreator->setNorm(norm);
    if(!m_surfaceCreator->execute()){
        m_surfaceCreator->sheet_cleaning();

        dualsheet.setSurface(m_surfaceCreator->getSurface());
        dual_sheets.push_back(dualsheet);

        for(auto const &t_v : m_surfaceCreator->getSurface()){
            (*m_block)[t_v.first] = -1;
            for(auto vec : t_v.second){
                colorTet(t_v.first,vec,sheet_ID);
            }
        }



        dualsheet.buildSurface(m_surfaceCreator->buildSurfaceSheet());

        //m_surfaceCreator->reset();
        //ASheetID = -1;
    }else{
        m_surfaceCreator->sheet_cleaning();
        dualsheet.setSurface(m_surfaceCreator->getSurface());
        dual_sheets.push_back(dualsheet);



        for(auto const &t_v : m_surfaceCreator->getSurface()){
            (*m_block)[t_v.first] = -1;
            for(auto vec : t_v.second){
                colorTet(t_v.first,vec,sheet_ID);
            }
        }



        dualsheet.buildSurface(m_surfaceCreator->buildSurfaceSheet());


        IGMeshIOService ioService_m(dualsheet.getSurfaceMesh());
        VTKWriter vtkWriter_m(&ioService_m);
        vtkWriter_m.setCellOptions(gmds::N|gmds::F);
        vtkWriter_m.setDataOptions(gmds::N|gmds::F);
        vtkWriter_m.write("/home/simon/Data/Results_debug/sheet_"+std::to_string(sheet_ID)+".vtk");

        sheet_ID++;
    }
    m_surfaceCreator->reset();
    return dualsheet.getSurface();
}
/*------------------------------------------------------------------------*/

std::vector<TCellID> DualBlockingSession::createBoundary(TCellID AID, int &ASheetID, Mesh* ASurface_mesh) {

    ASheetID = sheet_ID;

    DualSheet dualsheet(sheet_ID, ASurface_mesh);
    //cad::GeomMeshLinker linker_tmp = linkerG_T;

    std::vector<TCellID > nodes = m_mesh->get<Region>(AID).getIDs<Node>();
    for(auto n : nodes){
        if(linkerG_T.getGeomDim<Node>(n) == 3){
            m_boundaryCreator->setSurfaceID(linkerG_T.getGeomId<Node>(n));
        }
    }

    m_boundaryCreator->setSheetID(sheet_ID);
    m_boundaryCreator->setTetID(AID);
    if(!m_boundaryCreator->execute()){
        m_boundaryCreator->reset();
    }else {
        m_boundarySheets_to_surf.emplace(ASheetID, m_boundaryCreator->getSurfaceID());
        dualsheet.setSurface(m_boundaryCreator->getSurface());
        dualsheet.boundarySheet();
        dual_sheets.push_back(dualsheet);

        for(auto const &t_v : m_boundaryCreator->getSurface()){
            (*m_block)[t_v.first] = -1;
            for(auto vec:t_v.second){
                colorTet(t_v.first,vec,sheet_ID);
            }
        }
        dualsheet.buildSurface(m_boundaryCreator->buildSurfaceSheet());
        sheet_ID++;
    }

    for(auto const &t_v : dualsheet.getSurface()){
        (*m_block)[t_v] = -1;
    }
    m_boundaryCreator->reset();

    return dualsheet.getSurface();

}
/*------------------------------------------------------------------------*/

int DualBlockingSession::colorTet(TCellID AId, math::Vector3d AV, int ASheet) {

    double max  = 0;
    int axis=-1;
    math::Chart c = (*m_tetra_chart)[AId];

    for(int i = 0;i < 3; i++)
    {
        math::Vector3d vec = c[i];
        double product = vec.dot(AV);
        product = fabs(product);
        if(product > max)
        {
            max = product;
            axis = i;
        }
    }



    if(axis == 0){
        if((*m_sheet_X)[AId] == 0) {
            (*m_sheet_X)[AId] = ASheet;
        }
        else if((*m_sheet_X)[AId] == -1){

            std::cout<<"Singularity X in tet "<<AId<<std::endl;
            return 0;
        }
        else if(ASheet == 0){
            (*m_sheet_X)[AId] = 0;
        }
        else{
            std::cout<<"Sheet error in X "<<AId<<" with color "<<(*m_sheet_X)[AId]<<std::endl;
            return -1;
        }
    }
    else if(axis == 1){
        if((*m_sheet_Y)[AId] == 0) {
            (*m_sheet_Y)[AId] = ASheet;
        }
        else if((*m_sheet_Y)[AId] == -1){

            std::cout<<"Singularity Y in tet "<<AId<<std::endl;
            return 0;
        }
        else if(ASheet == 0){
            (*m_sheet_Y)[AId] = 0;
        }
        else{
            std::cout<<"Sheet error in Y "<<AId<<" with color "<<(*m_sheet_Y)[AId]<<std::endl;
            return -1;
        }
    }
    else {
        if((*m_sheet_Z)[AId] == 0) {
            (*m_sheet_Z)[AId] = ASheet;
        }
        else if((*m_sheet_Z)[AId] == -1){

            std::cout<<"Singularity Z in tet "<<AId<<std::endl;
            return 0;
        }
        else if(ASheet == 0){
            (*m_sheet_Z)[AId] = 0;
        }
        else{
            std::cout<<"Sheet error in Z "<<AId<<" with color "<<(*m_sheet_Z)[AId]<<std::endl;
            return -1;
        }

    }

    if(((*m_sheet_X)[AId] == (*m_sheet_Y)[AId] && (*m_sheet_X)[AId] !=0 && (*m_sheet_X)[AId] !=-1)||
       ((*m_sheet_Y)[AId] == (*m_sheet_Z)[AId] && (*m_sheet_Y)[AId] !=0 && (*m_sheet_Y)[AId] !=-1)||
       ((*m_sheet_X)[AId] == (*m_sheet_Z)[AId] && (*m_sheet_X)[AId] !=0 && (*m_sheet_X)[AId] !=-1)){

        std::cout<<"Tet "<<AId<<" in error : "<<std::endl;
        std::cout<<"\t X: "<<(*m_sheet_X)[AId]<<std::endl;
        std::cout<<"\t Y: "<<(*m_sheet_Y)[AId]<<std::endl;
        std::cout<<"\t Z: "<<(*m_sheet_Z)[AId]<<std::endl;
    }

    return 0;
}
/*------------------------------------------------------------------------*/
bool DualBlockingSession::isColoredIn(TCellID AId, int ASheet)
{
    return (*m_sheet_X)[AId] == ASheet || (*m_sheet_Y)[AId] == ASheet || (*m_sheet_Z)[AId] == ASheet;
}
/*------------------------------------------------------------------------*/
bool DualBlockingSession::isColored(TCellID AId)
{
    return (*m_sheet_X)[AId] != 0 || (*m_sheet_Y)[AId] != 0 || (*m_sheet_Z)[AId] != 0;
}
/*------------------------------------------------------------------------*/
int DualBlockingSession::getColor(gmds::TCellID AID, gmds::math::Vector3d AV) {

    double max  = 0;
    int axis=-1;
    math::Chart c = (*m_tetra_chart)[AID];

    for(int i = 0;i < 3; i++)
    {
        math::Vector3d vec = c[i];
        double product = vec.dot(AV);
        product = fabs(product);
        if(product > max)
        {
            max = product;
            axis = i;
        }
    }

    if(axis == 0){
        return (*m_sheet_X)[AID];
    }else if(axis == 1){
        return (*m_sheet_Y)[AID];
    }else{
        return (*m_sheet_Z)[AID];
    }
}

/*------------------------------------------------------------------------*/

bool DualBlockingSession::checkValidity() {

    /*bool dual_validity = !dual_sheet.empty() || !dual_bounder.empty();

    for(auto s : dual_sheet){
        std::map<TCellID,math::Vector3d> sheet_tets = s.getSurface();

        bool sheet_validity = false;

        for(auto tet : sheet_tets){
            if((*m_sheet_X)[tet.first] != 0 && (*m_sheet_Y)[tet.first] != 0 && (*m_sheet_Z)[tet.first] != 0  ){
                sheet_validity = true;
            }
        }

        if(!sheet_validity){
            std::cout<<"Sheet n°"<<s.getID()<<" is not valid, it's not intersected enough"<<std::endl;
            return false;
        }
    }

    for(auto b : dual_bounder){
        std::map<TCellID,math::Vector3d> sheet_tets = b.getSurface();

        bool sheet_validity = false;

        for(auto tet : sheet_tets){
            if((*m_sheet_X)[tet.first] != 0 && (*m_sheet_Y)[tet.first] != 0 && (*m_sheet_Z)[tet.first] != 0  ){
                sheet_validity = true;
            }
        }

        if(!sheet_validity){
            //std::cout<<"Sheet n°"<<b.getID()<<" is not valid, it's not intersected enough"<<std::endl;
            return false;
        }
    }

    return dual_validity;*/
}
/*------------------------------------------------------------------------*/
math::Chart DualBlockingSession::computeChart(std::vector<Node> ANodes, math::Point AP)
{

    math::Vector3d new_vectors[3];
    math::Chart charts[4];
    std::vector<math::Point> points;
    points.resize(4);
    for(int i = 0; i<ANodes.size(); i++){
        points[i] = ANodes[i].getPoint();
        charts[i] = (*m_vertex_chart)[ANodes[i].id()];
    }

    TCoord coords[4];
    math::Point::computeBarycentric(points[0], points[1], points[2], points[3], AP,
                                    coords[0], coords[1], coords[2],coords[3]);
    for(int i = 0; i<3; i++){
        math::Vector3d ref_i = charts[0][i];
        math::Vector3d vector = coords[0]*ref_i;

        for(int j = 1; j<4; j++){
            math::Vector3d v_j = closestComponentVector(ref_i,charts[j]);
            vector += coords[j]*v_j;
        }
        new_vectors[i] = vector;
    }
    math::Chart new_Chart(new_vectors[0],new_vectors[1],new_vectors[2]);
    return new_Chart;
}
/*----------------------------------------------------------------------------*/

std::vector<TCellID> DualBlockingSession::getSheet(int ASheet) {

    int i;

    for(i = 0; i<dual_sheets.size();i++){
        if(dual_sheets[i].getID() == ASheet){
            break;
        }
    }

    return dual_sheets[i].getSurface();
}
/*----------------------------------------------------------------------------*/
const math::Vector3d DualBlockingSession::closestComponentVector(const math::Vector3d& AV,
                                                                 math::Chart& AChart)
{
    //std::cout<<" ================= Finding closest component vector ================= "<<std::endl;

    math::Vector3d closest_vec;
    double max  = 0;
    for(int i = 0;i < 3; i++)
    {
        math::Vector3d vec = AChart[i];
        double product = vec.dot(AV);
        product = fabs(product);
        if(product > max)
        {
            max = product;
            closest_vec = vec;
        }
    }
    if( closest_vec.dot(AV) < 0 )
        closest_vec = -1*closest_vec;
    return closest_vec;
}
/*----------------------------------------------------------------------------*/

bool DualBlockingSession::colorDual(int &ANbZones) {

    m_dual_zones.clear();
    m_zone_border.clear();

    m_hmesh->clear();

    for(auto r : m_mesh->regions()){
        if((*m_block)[r] > 0){
            (*m_block)[r] = 0;
        }
    }
    for(auto r : m_mesh->regions()){
        if(m_part->value(r) > 1)
            m_block->set(r,-1);
    }

    std::cout<<"start dual zone"<<std::endl;

    std::clock_t start;
    double duration;

    int zone_id = 1;

    std::set<TCellID> zone;
    std::set<Node> nodes;

    std::set<TCellID> all_tet_in_zone;

    start = std::clock();
    std::set<std::pair<int,int>> dual_zone_classifications;

    for (auto t : m_mesh->regions()) {
        if(m_part->value(t) <= 1) {
            int type = 4;
            bool sing = false;
            TCellID sing_id = NullID;
            if ((*m_block)[t] == 0) {
                nodes.clear();
                std::set<int> borders;
                zone.clear();
                all_tet_in_zone.clear();
                //std::cout<<"first tet :"<<t<<std::endl;

                (*m_block)[t] = zone_id;

                std::vector<Node> tet_nodes = m_mesh->get<Region>(t).get<Node>();

                all_tet_in_zone.insert(t);

                for (auto n : tet_nodes) {

                    int type_tmp = 4;

                    if (linkerG_T.getGeomDim<Node>(n.id()) == cad::GeomMeshLinker::LINK_POINT) {
                        //std::cout << "link to point" << std::endl;
                        type_tmp = 1;
                        dual_zone_classifications.emplace(1,linkerG_T.getGeomId<Node>(n.id()));
                    } else if (linkerG_T.getGeomDim<Node>(n.id()) == cad::GeomMeshLinker::LINK_CURVE) {
                        //std::cout << "link to curve" << std::endl;
                        type_tmp = 2;
                        dual_zone_classifications.emplace(2,linkerG_T.getGeomId<Node>(n.id()));
                    } else if (linkerG_T.getGeomDim<Node>(n.id()) == cad::GeomMeshLinker::LINK_SURFACE) {
                        //std::cout << "link to surface" << std::endl;
                        type_tmp = 3;
                        dual_zone_classifications.emplace(3,linkerG_T.getGeomId<Node>(n.id()));
                    } else if (linkerG_T.getGeomDim<Node>(n.id()) == cad::GeomMeshLinker::LINK_VOLUME) {
                        //std::cout << "link to volume" << std::endl;
                        type_tmp = 4;
                        dual_zone_classifications.emplace(0,linkerG_T.getGeomId<Node>(n.id()));
                    } else {
                        //std::cout << "No link" << std::endl;
                    }


                    if (type_tmp < type && type_tmp != 0) {
                        //std::cout << "type changed from " << type << " to " << type_tmp
                        //          << std::endl;
                        //std::cout << "From tet " << t << std::endl;
                        type = type_tmp;
                        nodes.clear();

                        nodes.insert(n);
                    } else if (type_tmp == type) {
                        nodes.insert(n);
                    }


                    if ((*m_sing)[t] > 0) {
                        //std::cout << "sing found " << t << std::endl;
                        sing = true;
                        sing_id = t;
                    }
                }
                std::vector<Edge> edges = m_mesh->get<Region>(t).get<Edge>();

                //std::cout<<"edge size "<<m_mesh->get<Region>(t).get<Edge>().size()<<std::endl;

                for (auto e : edges) {
                    std::vector<TCellID> ids = e.getIDs<Region>();
                    for(auto id : ids){
                        if(m_part->value(id) <= 1)
                            zone.insert(ids.begin(), ids.end());

                    }
                }
                //std::cout<<"zone size "<<zone.size()<<std::endl;
                while (!zone.empty()) {
                    std::vector<TCellID> zone_vec;
                    zone_vec.insert(zone_vec.begin(), zone.begin(), zone.end());

                    zone.clear();

                    for (auto t_zone : zone_vec) {
                        if(m_part->value(t) <= 1) {
                            if ((*m_block)[t_zone] == 0) {
                                (*m_block)[t_zone] = zone_id;

                                all_tet_in_zone.insert(t_zone);

                                tet_nodes = m_mesh->get<Region>(t_zone).get<Node>();


                                for (auto n_zone : tet_nodes) {
                                    bool continue_cause_cut = false;
                                    for(auto t_cut : n_zone.getIDs<Region>()){
                                        if(m_part->value(t_cut) > 1){
                                            continue_cause_cut = true;
                                        }
                                    }
                                    if(continue_cause_cut) continue;

                                    int type_tmp = 4;

                                    if (linkerG_T.getGeomDim<Node>(n_zone.id()) == cad::GeomMeshLinker::LINK_POINT) {
                                        //std::cout << "link to point in tet "<<t_zone<<" node "<<n_zone.id() << std::endl;
                                        type_tmp = 1;
                                        dual_zone_classifications.emplace(1,linkerG_T.getGeomId<Node>(n_zone.id()));
                                    } else if (linkerG_T.getGeomDim<Node>(n_zone.id()) ==
                                               cad::GeomMeshLinker::LINK_CURVE) {
                                        //std::cout << "link to curve" << std::endl;
                                        type_tmp = 2;
                                        dual_zone_classifications.emplace(2,linkerG_T.getGeomId<Node>(n_zone.id()));
                                    } else if (linkerG_T.getGeomDim<Node>(n_zone.id()) ==
                                               cad::GeomMeshLinker::LINK_SURFACE) {
                                        //std::cout << "link to surface" << std::endl;
                                        type_tmp = 3;
                                        dual_zone_classifications.emplace(3,linkerG_T.getGeomId<Node>(n_zone.id()));
                                    } else if (linkerG_T.getGeomDim<Node>(n_zone.id()) ==
                                               cad::GeomMeshLinker::LINK_VOLUME) {
                                        //std::cout << "link to volume" << std::endl;
                                        type_tmp = 4;
                                        dual_zone_classifications.emplace(0,linkerG_T.getGeomId<Node>(n_zone.id()));
                                    } else {
                                        //std::cout << "No link" << std::endl;
                                    }

                                    if (type_tmp < type && type_tmp != 0) {
                                        //std::cout << "type changed from " << type << " to "
                                        //         << type_tmp << std::endl;
                                        //std::cout<<"From tet "<<t_zone<<std::endl;
                                        type = type_tmp;
                                        nodes.clear();

                                        nodes.insert(n_zone);
                                    } else if (type_tmp == type) {
                                        nodes.insert(n_zone);
                                    }
                                }

                                if ((*m_sing)[t_zone] > 0) {

                                    sing = true;
                                    sing_id = t_zone;
                                }

                                std::vector<Edge> edges_zone = m_mesh->get<Region>(t_zone).get<Edge>();

                                for (auto e : edges_zone) {
                                    std::vector<TCellID> ids = e.getIDs<Region>();
                                    for(auto id : ids){
                                        if(m_part->value(id) <= 1)
                                            zone.insert(ids.begin(), ids.end());

                                    }
                                }
                            } else if ((*m_block)[t_zone] == -1) {
                                int sheet_id = -1;
                                if ((*m_sheet_X)[t_zone] != 0 && (*m_sheet_X)[t_zone] != -1) {
                                    sheet_id = (*m_sheet_X)[t_zone];
                                    borders.insert(sheet_id);
                                }
                                if ((*m_sheet_Y)[t_zone] != 0 && (*m_sheet_Y)[t_zone] != -1) {
                                    sheet_id = (*m_sheet_Y)[t_zone];
                                    borders.insert(sheet_id);
                                }
                                if ((*m_sheet_Z)[t_zone] != 0 && (*m_sheet_Z)[t_zone] != -1) {
                                    sheet_id = (*m_sheet_Z)[t_zone];
                                    borders.insert(sheet_id);
                                }
                            }
                        }
                    }
                }


                if (all_tet_in_zone.size() == 1) {
                    std::cout << "Only one tet in the dual zone" << std::endl;
                    TCellID id = *all_tet_in_zone.begin();

                    (*m_block)[id] = -1;
                    continue;
                }

                math::Point vertex;

                Node block_node;

                std::vector<Node> nodes_vec;
                nodes_vec.insert(nodes_vec.begin(), nodes.begin(), nodes.end());
                //std::cout << "size " << nodes_vec.size() << std::endl;
                //std::cout << "type " << type << std::endl;

                if (type == 1) {
                    //Vertex zone
                    if (nodes_vec.size() != 1) {
                        //plusieurs sommets geometriques ou une singularite dans la zone
                        //ERROR
                        std::cout << "More than one geom point in the dual zone " << zone_id << std::endl;
                        //return false;
                    } else if (sing) {
                        //ERROR
                        std::cout << "More than one sing in the dual zone" << std::endl;
                        std::cout << sing_id << std::endl;
                        //return false;
                    }
                    vertex = nodes_vec[0].getPoint();

                    block_node = m_hmesh->newNode(vertex);
                    std::cout << "Node id in hmesh " << block_node.id() << std::endl;
                    linkerH_G.linkToPoint(block_node, linkerG_T.getGeomId<Node>(nodes_vec[0].id()));
                    //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                    //std::cout << "Corner point created at "<<block_node.getPoint() << std::endl;


                } else if (type == 2) {
                    //Edge zone
                    //if (sing) {

                    /*Region sing_tet = m_mesh->get<Region>(sing_id);
                    std::cout<<"sing_id : "<<sing_id<<std::endl;
                    vertex = sing_tet.center();

                    std::vector<cad::GeomSurface *> surfs;
                    manager.getSurfaces(surfs);
                    for (auto s : surfs) {
                        if (s->id() == linkerG_T.getGeomId<Node>(nodes_vec[0].id())) {
                            s->project(vertex);
                        }
                    }

                    block_node = m_hmesh->newNode(vertex);
                    linkerH_G.linkToSurface(block_node, linkerG_T.getGeomId<Node>(nodes_vec[0].id()));
                    m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                    std::cout << "Edge sing point created at " << block_node.getPoint() << std::endl;

                    //singularite dans la zone au bord
                    //ERROR
                    std::cout << "Sing in edge dual zone "<<zone_id <<" on tet "<<sing_id<< std::endl;
                    //return false;*/
                    //}else {

                    //std::cout << "Nb nodes " << nodes_vec.size() << std::endl;

                    vertex = nodes_vec[0].getPoint();
                    for (int i = 1; i < nodes_vec.size(); i++) {
                        vertex.setXYZ(vertex.X() + nodes_vec[i].getPoint().X(),
                                      vertex.Y() + nodes_vec[i].getPoint().Y(),
                                      vertex.Z() + nodes_vec[i].getPoint().Z());
                    }
                    vertex.setXYZ(vertex.X() / nodes_vec.size(),
                                  vertex.Y() / nodes_vec.size(),
                                  vertex.Z() / nodes_vec.size());

                    //std::cout << "point before projection " << vertex << std::endl;

                    std::vector<cad::GeomCurve *> curves;
                    manager.getCurves(curves);
                    for (auto c : curves) {
                        if (c->id() == linkerG_T.getGeomId<Node>(nodes_vec[0].id())) {
                            c->project(vertex);
                        }
                    }

                    //std::cout << "point after projection " << vertex << std::endl;

                    block_node = m_hmesh->newNode(vertex);
                    linkerH_G.linkToCurve(block_node, linkerG_T.getGeomId<Node>(nodes_vec[0].id()));
                    //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                    //std::cout << "Edge point created at " << block_node.getPoint() << std::endl;
                    //}

                } else if (type == 3) {
                    //Surface zone
                    //sing = false;
                    if (sing) {
                        //singularite dans la zone le point du bloc est construit dessus

                        Region sing_tet = m_mesh->get<Region>(sing_id);
                        //std::cout<<"sing_id : "<<sing_id<<std::endl;
                        vertex = sing_tet.center();

                        std::vector<cad::GeomSurface *> surfs;
                        manager.getSurfaces(surfs);
                        for (auto s : surfs) {
                            if (s->id() == linkerG_T.getGeomId<Node>(nodes_vec[0].id())) {
                                s->project(vertex);
                            }
                        }

                        block_node = m_hmesh->newNode(vertex);
                        linkerH_G.linkToSurface(block_node, linkerG_T.getGeomId<Node>(nodes_vec[0].id()));
                        //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                        //std::cout << "Face sing point created at " << block_node.getPoint() << std::endl;



                    } else {

                        vertex = nodes_vec[0].getPoint();
                        for (int i = 1; i < nodes_vec.size(); i++) {
                            vertex.setXYZ(vertex.X() + nodes_vec[i].getPoint().X(),
                                          vertex.Y() + nodes_vec[i].getPoint().Y(),
                                          vertex.Z() + nodes_vec[i].getPoint().Z());
                        }
                        vertex.setXYZ(vertex.X() / nodes_vec.size(),
                                      vertex.Y() / nodes_vec.size(),
                                      vertex.Z() / nodes_vec.size());

                        std::vector<cad::GeomSurface *> surfs;
                        manager.getSurfaces(surfs);
                        for (auto s : surfs) {
                            if (s->id() == linkerG_T.getGeomId<Node>(nodes_vec[0].id())) {
                                s->project(vertex);
                            }
                        }

                        block_node = m_hmesh->newNode(vertex);
                        linkerH_G.linkToSurface(block_node, linkerG_T.getGeomId<Node>(nodes_vec[0].id()));

                        //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                        //std::cout << "Face point created at " << block_node.getPoint() << std::endl;
                    }

                } else if (type == 4) {
                    //Volume zone
                    sing = false;
                    if (sing) {
                        //singularite dans la zone le point du bloc est construit dessus

                        Region sing_tet = m_mesh->get<Region>(sing_id);
                        vertex = sing_tet.center();
                        block_node = m_hmesh->newNode(vertex);
                        //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());

                        //std::cout << "Volume sing point created at " << block_node.getPoint() << std::endl;
                    } else {
                        vertex = math::Point(0, 0, 0);

                        for (auto v:nodes_vec) {
                            vertex = vertex + v.getPoint();
                        }
                        vertex = (1.0 / (double) nodes_vec.size()) * vertex;
                        block_node = m_hmesh->newNode(vertex);
                        //m_hmesh->newTriangle(block_node.id(), block_node.id(), block_node.id());
                        //std::cout << "Volume point created at " << block_node.getPoint() << std::endl;
                    }
                } else {
                    std::cout << "ERROR: zone " << zone_id << " of type " << type << std::endl;
                    exit(0);
                }

                bool good_zone = false;

                switch (type) {
                    case 1:
                        if (borders.size() == 3) {
                            good_zone = true;
                        } else {
                            //std::cout<<"ERROR: Point zone "<<zone_id<<" touch "<<borders.size()<<" dual sheets instead of 3"<<std::endl;
                            good_zone = true;
                        }
                        break;
                    case 2:
                        if (borders.size() == 4) {
                            good_zone = true;
                        } else {
                            //std::cout<<"ERROR: Curve zone "<<zone_id<<" touch "<<borders.size()<<" dual sheets instead of 4"<<std::endl;
                            good_zone = true;
                        }
                        break;
                    case 3:
                        if (borders.size() == 5) {
                            good_zone = true;
                        } else if (borders.size() == 3 && sing) {
                            good_zone = true;
                        } else if (borders.size() == 6 && sing) {
                            good_zone = true;
                        } else {
                            //std::cout<<"ERROR: Surface zone "<<zone_id<<" touch "<<borders.size()<<" dual sheets instead of 5"<<std::endl;
                            good_zone = true;
                        }
                        break;
                    case 4:
                        if (borders.size() == 6) {
                            good_zone = true;
                        } else if (borders.size() == 5 && sing) {
                            good_zone = true;
                        } else if (borders.size() == 7 && sing) {
                            good_zone = true;
                        } else {
                            //std::cout<<"ERROR: Volume zone "<<zone_id<<" touch "<<borders.size()<<" dual sheets instead of 6"<<std::endl;
                            good_zone = true;
                        }
                        break;
                    default:
                        break;
                }
                if (borders.size() <= 2) {
                    //good_zone = false;
                }
                if (good_zone) {
                    std::vector<int> borders_vec;
                    borders_vec.insert(borders_vec.begin(), borders.begin(), borders.end());
                    m_zone_border[zone_id] = borders_vec;

                    std::vector<std::pair<int,int>> zone_classification;
                    zone_classification.reserve(dual_zone_classifications.size());
                    dual_zone_classifications.clear();


                    if(m_mesh->isMarked<Region>(t,mark_ghost)){
                        //On est obligé d'avoir un traitement différent parce que la zone duale est définie sur des tet
                        // plats créés sur le bord
                        for (auto tet : all_tet_in_zone) {
                            Region r = m_mesh->get<Region>(tet);
                            // Pour chaque noeud de la zone duale on va récuperer les tet adjacents, ce qui simule le
                            // fait qu'on rajoute une couche de mailles entre le bord de la géométrie et le maillage
                            // ATTENTION : aucune preuve formelle qu'il n'y ait pas d'erreurs
                            for (auto n : r.get<Node>()) {
                                for(auto r_n : n.get<Region>()){
                                    for(auto n_r_n : r_n.get<Node>()){

                                        //Redondance de tet, méthode brutalle
                                        if (linkerG_T.getGeomDim(n_r_n) == 2) {
                                            dual_zone_classifications.emplace(2, linkerG_T.getGeomId(n_r_n));
                                        } else if (linkerG_T.getGeomDim(n_r_n) == 3) {
                                            dual_zone_classifications.emplace(3, linkerG_T.getGeomId(n_r_n));
                                        } else if (linkerG_T.getGeomDim(n_r_n) == 0) {
                                            dual_zone_classifications.emplace(0, linkerG_T.getGeomId(n_r_n));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else {
                        for (auto tet : all_tet_in_zone) {
                            Region r = m_mesh->get<Region>(tet);
                            for (auto n : r.get<Node>()) {
                                if (linkerG_T.getGeomDim(n) == 2) {
                                    dual_zone_classifications.emplace(2, linkerG_T.getGeomId(n));
                                } else if (linkerG_T.getGeomDim(n) == 3) {
                                    dual_zone_classifications.emplace(3, linkerG_T.getGeomId(n));
                                } else if (linkerG_T.getGeomDim(n) == 0) {
                                    dual_zone_classifications.emplace(0, linkerG_T.getGeomId(n));
                                }
                            }
                        }
                    }
                    std::cout<<"Zone id "<<zone_id<<std::endl;
                    for(auto geom_info : dual_zone_classifications){
                        zone_classification.push_back(geom_info);
                        std::cout<<"Dim "<<geom_info.first<<", ID "<<geom_info.second<<std::endl;
                    }
                    m_dual_zones_classifications[zone_id] = zone_classification;

                    m_dual_zones.emplace(zone_id, block_node);

                    zone_id++;
                }
            } else {
                continue;
            }
        }
    }

    duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;

    std::cout<<"nb point "<<m_hmesh->getNbNodes()<<std::endl;
    std::cout<<"nb zone "<<m_dual_zones.size()<<std::endl;



    int nb_blocks = 0;


    std::clock_t start2;
    double duration2;

    start2 = std::clock();
    //nb_blocks = createBlock();
    duration2 = (std::clock()-start2)/(double)CLOCKS_PER_SEC;

    for(auto d : durations){
        // std::cout<<"duration d = "<<d<<std::endl;
    }


    //std::cout<<"| dual validity time = "<<duration<<"| nb dual zones = "<<m_dual_zones.size()<<"| blocks time = "<<duration2<<"| nb blocks = "<<nb_blocks<<std::endl;

    ANbZones=zone_id-1;

    /*std::vector<int> bad_dual_region = checkDualRegionsValidity();
    if(!bad_dual_region.empty()){

        for(auto id : bad_dual_region){
            std::cout<<"bad dual region "<<id<<std::endl;
        }

        for(auto r : m_mesh->regions()){
            if(std::find(bad_dual_region.begin(),bad_dual_region.end(),m_block->value(r)) != bad_dual_region.end()){
                m_block->set(r,-2);
            }
        }
    }*/


    IGMeshIOService ioService_m(m_mesh);
    VTKWriter vtkWriter_m(&ioService_m);
    vtkWriter_m.setCellOptions(gmds::N|gmds::R);
    vtkWriter_m.setDataOptions(gmds::N|gmds::R);
    vtkWriter_m.write("/home/simon/Data/Results_debug/blocks.vtk");

    return checkDualRegionsValidity();
}
/*----------------------------------------------------------------------------*/

void DualBlockingSession::classifyZone(int AID) {

    for (auto t : m_mesh->regions()) {
        if((*m_sheet_X)[t] > 0 || (*m_sheet_Y)[t] > 0 || (*m_sheet_Z)[t] > 0){
            (*m_block)[t] = -1;
        }
    }
}
/*----------------------------------------------------------------------------*/
int DualBlockingSession::createBlock() {

    std::cout<<"start blocks"<<std::endl;

    std::map<int,std::vector<int>> block_dual;

    std::map<std::tuple<int,int,int,int,int,int,int,int>,int> blocks;

    int block_id = 0;

    //=========================================
    //HACK
    for(auto r : m_mesh->regions()){
        m_part->set(r, 1);
    }
    //=========================================

    for(auto r : m_mesh->regions()){
        if(m_part->value(r) == 1) {
            //get a tet who's the intersection of 3 sheet surface
            if ((*m_sheet_X)[r] > 0 && (*m_sheet_Y)[r] > 0 && (*m_sheet_Z)[r] > 0) {
                std::cout<<"Block"<<std::endl;


                int X = (*m_sheet_X)[r];
                int Y = (*m_sheet_Y)[r];
                int Z = (*m_sheet_Z)[r];

                std::set<int> colors;

                std::vector<int> dual;

                dual.push_back(X);
                dual.push_back(Y);
                dual.push_back(Z);

                for (auto z_b : m_zone_border) {

                    if (std::find(z_b.second.begin(), z_b.second.end(), X) != z_b.second.end()) {
                        if (std::find(z_b.second.begin(), z_b.second.end(), Y) != z_b.second.end()) {
                            if (std::find(z_b.second.begin(), z_b.second.end(), Z) != z_b.second.end()) {

                                colors.insert(z_b.first);
                            }
                        }
                    }
                }
                if (colors.size() != 8) {
                    std::cout << "Tet " << r << " for block " << block_id << std::endl;
                    std::cout << "Nb zones from triplet " << colors.size() << std::endl;
                }
                std::set<int> colors_dist;
                if (colors.size() > 8) {

                    std::vector<Node> nodes = m_mesh->get<Region>(r).get<Node>();
                    std::set<Region> surrounding;
                    for (auto n : nodes) {
                        std::vector<Region> tmp = n.get<Region>();
                        surrounding.insert(tmp.begin(), tmp.end());
                    }
                    std::vector<Region> surrounding_vec;
                    do {

                        surrounding_vec.insert(surrounding_vec.begin(), surrounding.begin(), surrounding.end());
                        surrounding.clear();

                        for (auto t : surrounding_vec) {
                            if(m_part->value(t.id()) == 1) {
                                if ((*m_block)[t.id()] != -1) {
                                    if (std::find(colors.begin(), colors.end(), (*m_block)[t.id()]) != colors.end()) {
                                        colors_dist.insert((*m_block)[t.id()]);
                                    }
                                } else {
                                    nodes = t.get<Node>();
                                    for (auto n : nodes) {
                                        std::vector<Region> tmp = n.get<Region>();
                                        surrounding.insert(tmp.begin(), tmp.end());
                                    }
                                }
                            }
                        }

                        surrounding_vec.clear();
                    } while (colors_dist.size() < 8);
                    colors = colors_dist;
                }
                if (colors.size() == 8) {
                    bool wrong_block = false;
                    std::vector<int> block;

                    for (auto c : colors) {
                        block.push_back(c);

                    }

                    if (!wrong_block) {
                        int size_b = blocks.size();

                        if (blocks.emplace(
                                std::tuple<int, int, int, int, int, int, int, int>(block[0], block[1], block[2],
                                                                                   block[3],
                                                                                   block[4], block[5], block[6],
                                                                                   block[7]),
                                block_id).second) {
                            block_dual.emplace(block_id, dual);
                            block_id++;
                        }
                    }
                }
            }

            //if(block_id == 1)
            //break;
        }
    }

    std::cout<<"NB bloc found "<<blocks.size()<<std::endl;

    int b_index = 0;
    for(auto b : blocks){


        int cpt_block = 0;
        std::map<int,std::vector<int>> corners;


        int id0 = -1;
        int id1 = -1;
        int id2 = -1;
        int id3 = -1;
        int id4 = -1;
        int id5 = -1;
        int id6 = -1;
        int id7 = -1;


        bool block_good;
        bool to_much_it = false;

        do {
            cpt_block++;


            block_good = true;

            m_mesh->unmarkAll<Region>(mark_visited);


            std::vector<int> corner = createCorner(std::get<0>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<0>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<0>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<0>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<0>(b.first), b.first, block_dual[b.second]);
            }

            //std::cout << "Corner from zone " << std::get<0>(b.first) << ":" << std::endl;
            //std::cout<<"corner size "<<corner.size()<<std::endl;
            //std::cout << "\t" << corner[0] << std::endl;
            id1 = corner[0];
            //std::cout << "\t" << corner[1] << std::endl;
            id3 = corner[1];
            //std::cout << "\t" << corner[2] << std::endl;
            id4 = corner[2];

            corners.emplace(std::get<0>(b.first), corner);


            corner = createCorner(std::get<1>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<1>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<1>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<1>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<1>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<1>(b.first), corner);


            corner = createCorner(std::get<2>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<2>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<2>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<2>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<2>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<2>(b.first), corner);


            corner = createCorner(std::get<3>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<3>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<3>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<3>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<3>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<3>(b.first), corner);


            corner = createCorner(std::get<4>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<4>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<4>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<4>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<4>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<4>(b.first), corner);



            corner = createCorner(std::get<5>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<5>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<5>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<5>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<5>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<5>(b.first), corner);


            corner = createCorner(std::get<6>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<6>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<6>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<6>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<6>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<6>(b.first), corner);


            corner = createCorner(std::get<7>(b.first), b.first, block_dual[b.second]);
            //std::cout<<"Corner size in createBlock "<<corner.size()<<std::endl;
            if(corner.size() < 3){
                //std::cout<<"corner wrong"<<std::endl;
                block_good = false;

                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<7>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }

                corner = createCorner(std::get<7>(b.first), b.first, block_dual[b.second]);
                //break;
            }else if(corner.size() == 4){
                for(auto r : m_mesh->regions()){
                    if((*m_block)[r] == std::get<7>(b.first)){
                        m_mesh->unmark<Region>(r,mark_visited);
                    }
                }
                corner = createCorner(std::get<7>(b.first), b.first, block_dual[b.second]);
            }

            corners.emplace(std::get<7>(b.first), corner);
            //}

            id0 = std::get<0>(b.first);


            std::vector<int> corner1 = corners[corners[id0][0]];
            std::vector<int> corner3 = corners[corners[id0][1]];


            for (auto i1 : corner1) {
                //std::cout << "i1 = " << i1 << std::endl;
                for (auto i3 : corner3) {
                    //std::cout << "\ti3 = " << i3 << std::endl;
                    if (i1 != std::get<0>(b.first) && i1 == i3) {
                        id2 = i1;
                    }
                }
            }
            for (auto i1 : corner1) {
                if (i1 != id0 && i1 != id2) {
                    id5 = i1;
                }
            }
            for (auto i2 : corners[id2]) {
                if (i2 != id1 && i2 != id3) {
                    id6 = i2;
                }
            }
            for (auto i3 : corner3) {
                if (i3 != id0 && i3 != id2) {
                    id7 = i3;
                }
            }

            block_good = !(id2 == -1 || id5 == -1 || id6 == -1 || id7 == -1);

            std::vector<int> ids;
            ids.resize(8);
            ids[0]=id0;
            ids[1]=id1;
            ids[2]=id2;
            ids[3]=id3;
            ids[4]=id4;
            ids[5]=id5;
            ids[6]=id6;
            ids[7]=id7;

            bool zone_good;

            int to_reset = 0;

            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    if (j == i) continue;
                    if (ids[j] == ids[i]) {
                        //exit(1);
                        std::cout << "ids[" << j << "] = " << ids[j] << " == ids[" << i << "] = " << ids[i]
                                  << std::endl;
                        std::cout << "hex = " << id0 << "," << id1 << "," << id2 << "," << id3 << "," << id4 << ","
                                  << id5 << "," << id6 << "," << id7 << std::endl;
                        std::cout<<std::get<0>(b.first)<<","
                                 <<std::get<1>(b.first)<<","
                                 <<std::get<2>(b.first)<<","
                                 <<std::get<3>(b.first)<<","
                                 <<std::get<4>(b.first)<<","
                                 <<std::get<5>(b.first)<<","
                                 <<std::get<6>(b.first)<<","
                                 <<std::get<7>(b.first)<<","<<std::endl;
                        exit(1);
                        block_good = false;
                    }
                }
            }


        }while(!block_good);



        m_mesh->unmarkAll<Region>(mark_visited);

        m_hmesh->newHex(id0-1,id1-1,id2-1,id3-1,id4-1,id5-1,id6-1,id7-1);
        b_index++;

        //break;
    }


    Variable<int>* n_on_point = m_hmesh->getOrCreateVariable<int,GMDS_NODE>("n_on_point");
    Variable<int>* n_on_curve = m_hmesh->getOrCreateVariable<int,GMDS_NODE>("n_on_curve");
    Variable<int>* n_on_surface = m_hmesh->getOrCreateVariable<int,GMDS_NODE>("n_on_surface");
    for(auto n : m_hmesh->nodes()){
        n_on_point->set(n,-1);
        n_on_curve->set(n,-1);
        n_on_surface->set(n,-1);
        if(linkerH_G.getGeomDim<Node>(n) == 1){
            n_on_point->set(n, linkerH_G.getGeomId<Node>(n));
        }else if(linkerH_G.getGeomDim<Node>(n) == 2){
            n_on_curve->set(n, linkerH_G.getGeomId<Node>(n));
        }else if(linkerH_G.getGeomDim<Node>(n) == 3){
            n_on_surface->set(n, linkerH_G.getGeomId<Node>(n));
        }
    }


    MeshDoctor doch(m_hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    //doch.buildN2E(m_hmesh->getModel());
    doch.updateUpwardConnectivity();


    IGMeshIOService ioService(m_hmesh);
    VTKWriter vtkWriter(&ioService);
    //vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
    //vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::R);
    //vtkWriter.write("/home/simon/Data/Results_debug/hexa.vtk");

    classifyEdges(1);
    classifyBlockFaces();

    //correctionBlockEdgesClassification();

    vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::R);
    vtkWriter.write("/home/simon/Data/Results_debug/hexa.vtk");

    //smoothBlocks();

    /*vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write("/home/simon/Data/Results_debug/hexa_smooth.vtk");*/

    interpolationV2(10);

    //interpolationV2(5);

    return blocks.size();
}

/*----------------------------------------------------------------------------*/

void DualBlockingSession::checkContact() {

    /*for(auto s : dual_sheets){
        std::map<TCellID,math::Vector3d> sheet_tets = s.getSurface();

        for(auto tet : sheet_tets){

            Region t_sheet = m_mesh->get<Region>(tet.first);
            std::vector<Face> faces = t_sheet.get<Face>();
            for(auto f : faces){
                std::vector<Region> tets = f.get<Region>();
                for(auto t : tets){
                    int color = getColor(t.id(),tet.second);
                    if(color != 0 && color != s.getID()){

                        //Here the tetra tet is in contact with a tetra t from another sheet
                        math::Point center = t_sheet.center();
                        Node n = m_mesh->newNode(center);
                        for (auto f_sub : faces) {

                        }
                    }
                }
            }
        }
    }*/
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::smoothBlocks(){

    for(auto n : m_hmesh->nodes()){
        Node node = m_hmesh->get<Node>(n);
        std::cout<<"Node dim = "<<linkerH_G.getGeomDim<Node>(n)<<" id "<<linkerH_G.getGeomId<Node>(n)<<std::endl;
        std::cout<<"Nb edges for the node "<<node.nbEdges()<<std::endl;
        std::cout<<"Nb faces for the node "<<node.nbFAces()<<std::endl;
        std::cout<<"Nb regions for the node"<<node.nbRegions()<<std::endl;

    }

    /*for(int i = 0; i<10; i++) {
        for (auto n_id : m_hmesh->nodes()) {
            if (linkerH_G.getGeomDim<Node>(n_id) != 1) {

                Node n = m_hmesh->get<Node>(n_id);
                std::vector<math::Point> points;
                for (auto e : n.get<Edge>()) {
                    if (e.getIDs<Node>()[0] == n_id) {
                        points.push_back(e.get<Node>()[1].getPoint());
                    } else {
                        points.push_back(e.get<Node>()[0].getPoint());
                    }
                }
                math::Point p_n = math::Point::massCenter(points);

                if (linkerH_G.getGeomDim<Node>(n_id) == 2) {
                    cad::GeomCurve *c = manager.getCurve(linkerH_G.getGeomId<Node>(n_id));
                    c->project(p_n);
                } else if (linkerH_G.getGeomDim<Node>(n_id) == 3) {
                    cad::GeomSurface *s = manager.getSurface(linkerH_G.getGeomId<Node>(n_id));
                    s->project(p_n);
                }

                n.setPoint(p_n);
            }
        }
    }*/

    cad::GeomSmoother smoother(&linkerH_G);

    //smoother.smoothCurves(10);
    std::cout<<"curves"<<std::endl;
    //smoother.smoothSurfaces(10);
    std::cout<<"surfaces"<<std::endl;
    //smoother.smoothVolumes(10);
    std::cout<<"volumes"<<std::endl;

}
/*----------------------------------------------------------------------------*/

void DualBlockingSession::interpolationV2(int Ait){
    TCoord min[3],max[3];
    std::map<TCellID ,int> edges_sub;
    //manager.getVolume(1)->computeBoundingBox(min,max);
    //double distance = sqrt((min[0]-max[0])*(min[0]-max[0])+(min[1]-max[1])*(min[1]-max[1])+(min[2]-max[2])*(min[2]-max[2]));

    /*double min_length_hexa = 999999999;
    for(auto e : m_hmesh->edges()){
        Edge edge = m_hmesh->get<Edge>(e);
        if(edge.length() < min_length_hexa) min_length_hexa = edge.length();
    }

    double length = distance/30.0;
    for(auto e : m_hmesh->edges()){
        Edge edge = m_hmesh->get<Edge>(e);
        int nb_sub = (edge.length()/length)>=1?lround(edge.length())/length:1;
        std::cout<<"nb sub "<<nb_sub<<std::endl;
        edges_sub[e] = nb_sub;
    }

    bool changement = true;
    while(changement) {
        changement = false;
        for (auto f : m_hmesh->faces()) {
            Face face = m_hmesh->get<Face>(f);
            TCellID e0 = face.getIDs<Edge>()[0];
            TCellID e1 = face.getIDs<Edge>()[1];
            TCellID e2 = face.getIDs<Edge>()[2];
            TCellID e3 = face.getIDs<Edge>()[3];

            int indice_opp = 0;
            for(int i = 1; i<4; i++){
                if(face.get<Edge>()[i].getIDs<Node>()[0] != face.get<Edge>()[0].getIDs<Node>()[0] && face.get<Edge>()[i].getIDs<Node>()[0] != face.get<Edge>()[0].getIDs<Node>()[1] &&
                        face.get<Edge>()[i].getIDs<Node>()[1] != face.get<Edge>()[0].getIDs<Node>()[0] && face.get<Edge>()[i].getIDs<Node>()[1] != face.get<Edge>()[0].getIDs<Node>()[1]){
                    indice_opp = i;
                }
            }
            if(indice_opp == 1){
                e2 = face.getIDs<Edge>()[1];

                e1 = face.getIDs<Edge>()[2];
                e3 = face.getIDs<Edge>()[3];
            }else if(indice_opp == 2){
                e2 = face.getIDs<Edge>()[2];

                e1 = face.getIDs<Edge>()[1];
                e3 = face.getIDs<Edge>()[3];
            }else if(indice_opp == 3){
                e2 = face.getIDs<Edge>()[3];

                e1 = face.getIDs<Edge>()[1];
                e3 = face.getIDs<Edge>()[2];
            }

            if (edges_sub[e0] < edges_sub[e2]) {
                changement = true;
                edges_sub[e2] = edges_sub[e0];
            } else if(edges_sub[e0] > edges_sub[e2]){
                changement = true;
                edges_sub[e0] = edges_sub[e2];
            }

            if (edges_sub[e1] < edges_sub[e3]) {
                changement = true;
                edges_sub[e3] = edges_sub[e1];
            } else if(edges_sub[e1] > edges_sub[e3]){
                changement = true;
                edges_sub[e1] = edges_sub[e3];
            }
        }
    }


    for(auto e : m_hmesh->edges()){
        std::cout<<"nb sub "<<edges_sub[e]<<std::endl;
    }*/


    EdgeDiscrAlgo mapping_algo(m_hmesh);
    mapping_algo.init();
    edges_sub = mapping_algo.getEdgeDiscretization();
    for(auto e : m_hmesh->edges()){
        std::cout<<"edge "<<e<<" sub = "<<edges_sub[e]<<std::endl;
    }

    //manager.initFrom3DMesh(m_mesh);

    //linkerH_G.setGeometry(&manager);

    cad::GeomMeshLinker linker_hexa;
    linker_hexa.setGeometry(&manager);


    Mesh t_mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                          R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&t_mesh);
    MeshDoctor tdoc(&t_mesh);
    tdoc.buildFacesAndR2F();
    tdoc.buildEdgesAndX2E();
    tdoc.updateUpwardConnectivity();
    linker_hexa.setMesh(&t_mesh);

    double iteration = Ait;

    for(auto n : m_hmesh->nodes()){
        Node node = t_mesh.newNode(m_hmesh->get<Node>(n).getPoint());
        if(linkerH_G.getGeomDim<Node>(n) == 1){
            linker_hexa.linkToPoint(node,linkerH_G.getGeomId<Node>(n));
        } else if (linkerH_G.getGeomDim<Node>(n) == 2){
            linker_hexa.linkToCurve(node,linkerH_G.getGeomId<Node>(n));
        } else if (linkerH_G.getGeomDim<Node>(n) == 3){
            linker_hexa.linkToSurface(node,linkerH_G.getGeomId<Node>(n));
        }
    }


    std::map<TCellID,std::vector<int>> edge_discr;
    std::map<TCellID,std::map<int,std::vector<int>>> face_discr;

    std::vector<Region> blocks;

    m_hmesh->getAll(blocks);

    for(auto e : m_hmesh->edges()){
        //for(auto e : blocks[3].getIDs<Edge>()){
        std::cout<<"Edge discretization "<<e<<std::endl;

        Edge edge = m_hmesh->get<Edge>(e);

        Node first = edge.get<Node>()[0];
        Node second = edge.get<Node>()[1];

        int geom_dim_e1,geom_dim_e2;
        int geom_id_e1,geom_id_e2;

        bool border_edge = false;
        if((linkerH_G.getGeomDim(first) == 2 || linkerH_G.getGeomDim(first) == 3) &&
           ((linkerH_G.getGeomDim(second) == 2 || linkerH_G.getGeomDim(second) == 3))) {

            for (auto f : edge.get<Face>()) {
                if (f.get<Region>().size() == 1) {
                    border_edge = true;
                    std::cout<<"Border edge"<<std::endl;
                }
            }
        }

        int geom_id = -1;
        int type = 0;

        //Classification of the block edge
        /*if(border_edge){

            geom_id_e1 = linkerH_G.getGeomId(first);
            geom_id_e2 = linkerH_G.getGeomId(second);

            if(linkerH_G.getGeomDim(first) == 2 && linkerH_G.getGeomDim(second) == 2){
                if(geom_id_e1 != geom_id_e2) {
                    //Edge type = face
                    int gc1 = geom_id_e1;
                    int gc2 = geom_id_e2;

                    std::vector<cad::GeomCurve *> curves;
                    manager.getCurves(curves);

                    cad::GeomCurve *GC1;
                    cad::GeomCurve *GC2;

                    for (auto c : curves) {
                        if (c->id() == gc1) {
                            GC1 = c;
                        }
                        if (c->id() == gc2) {
                            GC2 = c;
                        }
                    }
                    geom_id = manager.getCommonSurface(GC1, GC2);
                    type = 3;
                    linkerH_G.linkEdgeToSurface(e,geom_id);
                }else{
                    //Edge type = edge
                    geom_id = geom_id_e1;
                    type = 2;
                    linkerH_G.linkEdgeToCurve(e,geom_id);
                }
            }else if(linkerH_G.getGeomDim(first) == 3 && linkerH_G.getGeomDim(second) == 3){
                if(geom_id_e1 != geom_id_e2) {
                    //Edge type = volume
                    type = 4;
                } else{
                    //Edge type = face
                    geom_id = geom_id_e1;
                    type = 3;
                    linkerH_G.linkEdgeToSurface(e,geom_id);
                }
            }else if(linkerH_G.getGeomDim(first) == 3){
                //Edge type = face
                geom_id = geom_id_e1;
                type = 3;
                linkerH_G.linkEdgeToSurface(e,geom_id);
            }else if(linkerH_G.getGeomDim(second) == 3){
                //Edge type = face
                geom_id = geom_id_e2;
                type = 3;
                linkerH_G.linkEdgeToSurface(e,geom_id);
            }

        }
        else {
            if (linkerH_G.getGeomDim(first) <= linkerH_G.getGeomDim(second)) {
                geom_dim_e1 = linkerH_G.getGeomDim(first);
                geom_dim_e2 = linkerH_G.getGeomDim(second);

                geom_id_e1 = linkerH_G.getGeomId(first);
                geom_id_e2 = linkerH_G.getGeomId(second);
            } else {
                geom_dim_e1 = linkerH_G.getGeomDim(second);
                geom_dim_e2 = linkerH_G.getGeomDim(first);

                geom_id_e1 = linkerH_G.getGeomId(second);
                geom_id_e2 = linkerH_G.getGeomId(first);
            }

            if (geom_dim_e1 == 1) {
                if (geom_dim_e2 == 1) {
                    //Edge type = edge
                    std::cout<<"Edge between two geom point"<<std::endl;
                    int gp1 = geom_id_e1;
                    int gp2 = geom_id_e2;

                    std::vector<cad::GeomPoint *> points;
                    manager.getPoints(points);

                    cad::GeomPoint *GP1;
                    cad::GeomPoint *GP2;

                    for (auto p : points) {
                        if (p->id() == gp1) {
                            GP1 = p;
                        }
                        if (p->id() == gp2) {
                            GP2 = p;
                        }
                    }

                    std::cout<<GP1->point()<<std::endl;
                    std::cout<<GP2->point()<<std::endl;

                    std::cout<<"NB Geom edge "<<manager.getNbCurves()<<std::endl;


                    geom_id = manager.getCommonCurve(GP1, GP2);
                    std::cout<<"Geom id = "<<geom_id<<std::endl;
                    type = 2;
                    linkerH_G.linkEdgeToCurve(e,geom_id);
                } else if (geom_dim_e2 == 2) {

                    int gp1 = geom_id_e1;
                    std::vector<cad::GeomPoint *> points;
                    manager.getPoints(points);

                    cad::GeomPoint *GP1;

                    for (auto p : points) {
                        if (p->id() == gp1) {
                            GP1 = p;
                        }
                    }

                    math::Point point_gp1 = GP1->point();

                    std::vector<cad::GeomCurve *> curves;
                    manager.getCurves(curves);

                    cad::GeomCurve *GC1;

                    for (auto c : curves) {
                        if (c->id() == geom_id_e2) {
                            GC1 = c;
                        }
                    }

                    math::Point closest = GC1->closestPoint(point_gp1);

                    std::cout << "point_gp1 = " << point_gp1 << std::endl;
                    std::cout << "closest = " << closest << std::endl;

                    if (closest.X() != point_gp1.X() || closest.Y() != point_gp1.Y() ||
                        closest.Z() != point_gp1.Z()) {
                        //Edge type = volume
                        std::cout << "point---edge" << std::endl;
                        type = 4;
                    } else {
                        //Edge type = edge
                        geom_id = geom_id_e2;
                        type = 2;
                        linkerH_G.linkEdgeToCurve(e,geom_id);
                    }

                } else if (geom_dim_e2 == 3) {
                    //Edge type = face
                    geom_id = geom_id_e2;
                    type = 3;
                    linkerH_G.linkEdgeToSurface(e,geom_id);
                } else {
                    //Edge type = volume
                    type = 4;
                }
            } else if (geom_dim_e1 == 2) {
                if (geom_dim_e2 == 2) {
                    if (geom_id_e1 == geom_id_e2) {
                        //Edge type = edge
                        geom_id = geom_id_e1;
                        type = 2;
                        linkerH_G.linkEdgeToCurve(e,geom_id);
                    } else {
                        int gc1 = geom_id_e1;
                        int gc2 = geom_id_e2;

                        std::vector<cad::GeomCurve *> curves;
                        manager.getCurves(curves);

                        cad::GeomCurve *GC1;
                        cad::GeomCurve *GC2;

                        for (auto c : curves) {
                            if (c->id() == gc1) {
                                GC1 = c;
                            }
                            if (c->id() == gc2) {
                                GC2 = c;
                            }
                        }

                        geom_id = manager.getCommonSurface(GC1, GC2);

                        if(geom_id==-1){
                            //Edge type = volume
                            type = 4;
                        }
                        else{
                            //Edge type = face
                            type = 3;
                            linkerH_G.linkEdgeToSurface(e,geom_id);
                        }
                    }
                } else if (geom_dim_e2 == 3) {
                    //Edge type = volume
                    type = 4;
                } else {
                    //Edge type = volume
                    type = 4;
                }
            } else if (geom_dim_e1 == 3) {
                if (geom_dim_e2 == 3) {
                    if (geom_id_e1 == geom_id_e2) {
                        //Edge type = face
                        geom_id = geom_id_e1;
                        type = 3;
                        linkerH_G.linkEdgeToSurface(e,geom_id);
                    }
                } else {
                    //Edge type = volume
                    type = 4;
                }
            }
        }*/

        math::Point p_first;
        math::Point p_second;

        std::vector<int> vec_allocator;
        vec_allocator.resize(edges_sub[e]+1);
        edge_discr[e] = vec_allocator;
        type = linkerH_G.getGeomDim<Edge>(e);
        geom_id = linkerH_G.getGeomId<Edge>(e);


        if(first.id() < second.id()){
            p_first = first.getPoint();
            p_second = second.getPoint();

            edge_discr[e][0] = first.id();
            edge_discr[e][edges_sub[e]] = second.id();
        }else{
            p_first = second.getPoint();
            p_second = first.getPoint();
            edge_discr[e][0] = second.id();
            edge_discr[e][edges_sub[e]] = first.id();
        }

        std::cout<<"First point "<<first.id()<<std::endl;
        std::cout<<"Second point "<<second.id()<<std::endl;

//t_mesh.newNode(p_first);
//t_mesh.newNode(p_second);

        math::Vector3d vector(p_first, p_second);
        math::Vector3d mini_vec = vector / 10;

        for (int i = 1; i < edges_sub[e]; i++) {
            //math::Point point(p1+mini_vec*i);
            math::Point point;
            point = (p_second * (double(i) / edges_sub[e]) + p_first * (1 - (double(i) / edges_sub[e])));
            //std::cout<<"Point: "<<point<<std::endl;
            //std::cout<<"point:"<<point<<std::endl;
            //std::cout<<"vector Y: "<<vector.Y()<<std::endl;
            if (type == 2) {
                std::vector<cad::GeomCurve *> curves;
                manager.getCurves(curves);
                for (auto c : curves) {
                    if (c->id() == geom_id) {
                        //std::cout<<"Geom id for projection "<<c->id()<<std::endl;
                        c->project(point);
                    }
                }
            } else if (type == 3) {
                std::vector<cad::GeomSurface *> surfs;
                manager.getSurfaces(surfs);
                for (auto s : surfs) {
                    if (s->id() == geom_id) {
                        //s->project(point);
                    }
                }
            }

            //Node inter_node = m_hmesh->newNode(point);

            Node inter_node = t_mesh.newNode(point);

            // t_mesh.newTriangle(inter_node,inter_node,inter_node);
            if (type == 2) {
                linker_hexa.linkToCurve(inter_node, geom_id);
            } else if (type == 3) {
                linker_hexa.linkToSurface(inter_node, geom_id);
            }

            edge_discr[e][i] = inter_node.id();
        }

        for(int i_p = 0; i_p<edges_sub[e]; i_p++){
            std::cout<<edge_discr[e][i_p]<<std::endl;
            //t_mesh.newTriangle(edge_discr[e][i_p],edge_discr[e][i_p],edge_discr[e][i_p+1]);

        }

        std::cout<<"nb points "<<t_mesh.getNbNodes()<<std::endl;
        std::cout<<"nb faces "<<t_mesh.getNbFaces()<<std::endl;

        //break;
    }

    for(auto f : m_hmesh->faces()) {

        Face face = m_hmesh->get<Face>(f);
        std::vector<Node> nodes = face.get<Node>();
        std::vector<Edge> edges = face.get<Edge>();

        int ex,ey;
        for(auto e : edges){
            std::cout<<"edge of face sub "<<edges_sub[e.id()]<<std::endl;
            if((e.getIDs<Node>()[0] == nodes[0].id() && e.getIDs<Node>()[1] == nodes[1].id()) ||
               (e.getIDs<Node>()[1] == nodes[0].id() && e.getIDs<Node>()[0] == nodes[1].id())){
                ex = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[1].id() && e.getIDs<Node>()[1] == nodes[2].id()) ||
                     (e.getIDs<Node>()[1] == nodes[1].id() && e.getIDs<Node>()[0] == nodes[2].id())){
                ey = e.id();
            }
        }
        int subx = edges_sub[ex];
        int suby = edges_sub[ey];


        std::vector<int> vec_allocator;
        vec_allocator.resize(edges_sub[ey]+1,0);
        for(int i = 0; i<=edges_sub[ex]; i++){
            face_discr[f][i] = vec_allocator;
        }/*
        for (auto tab : face_discr[f]) {
            std::vector<int> vec_allocator;
            vec_allocator.resize(edges_sub[ey]+1,0);
            tab.second = vec_allocator;
        }*/
        for(int i = 1; i < edges_sub[ex]; i++){
            for(int j = 1; j < edges_sub[ey]; j++){
                face_discr[f][i][j] = -1;
            }
        }

        std::cout << "Face = ";
        for (auto n : nodes) {
            std::cout << n.id() << " ";
        }
        std::cout << std::endl;

        for (int i = 0; i < 3; i++) {
            std::cout << "Nodes : " << nodes[i].id() << "," << nodes[i + 1].id() << std::endl;
            for (auto e: edges) {
                std::cout << "Edge " << edge_discr[e.id()][0] << "," << edge_discr[e.id()][edges_sub[e.id()]] << std::endl;
                if (nodes[i].id() == edge_discr[e.id()][0] &&
                    nodes[i + 1].id() == edge_discr[e.id()][edges_sub[e.id()]]) {
                    // + direction
                    std::cout << "+ direction" << std::endl;

                    switch (i) {
                        case 0:
                            for (int x = 0; x <= edges_sub[e.id()]; x++) {
                                face_discr[f][x].reserve(edges_sub[e.id()] + 1);
                                face_discr[f][x][0] = edge_discr[e.id()][x];
                                //std::cout << edge_discr[e.id()][x] << std::endl;
                            }
                            break;
                        case 1:
                            for (int y = 0; y <= edges_sub[e.id()]; y++) {
                                face_discr[f][edges_sub[ex]][y] = edge_discr[e.id()][y];
                                //std::cout << edge_discr[e.id()][y] << std::endl;
                            }
                            break;
                        case 2:
                            for (int x = 0; x <= edges_sub[e.id()]; x++) {
                                face_discr[f][x][edges_sub[ey]] = edge_discr[e.id()][edges_sub[e.id()] - x];
                                //std::cout << edge_discr[e.id()][edges_sub[e.id()] - x] << std::endl;
                            }
                            break;
                        default:
                            break;
                    }

                } else if (nodes[i].id() == edge_discr[e.id()][edges_sub[e.id()]] &&
                           nodes[i + 1].id() == edge_discr[e.id()][0]) {
                    // - direction
                    std::cout << "- direction" << std::endl;

                    switch (i) {
                        case 0:
                            for (int x = 0; x <= edges_sub[e.id()]; x++) {
                                face_discr[f][x].reserve(edges_sub[e.id()] + 1);
                                face_discr[f][x][0] = edge_discr[e.id()][edges_sub[e.id()] - x];
                                //std::cout << edge_discr[e.id()][edges_sub[e.id()] - x] << std::endl;
                            }
                            break;
                        case 1:
                            for (int y = 0; y <= edges_sub[e.id()]; y++) {
                                face_discr[f][edges_sub[ex]][y] = edge_discr[e.id()][edges_sub[e.id()] - y];
                                //std::cout << edge_discr[e.id()][edges_sub[e.id()] - y] << std::endl;
                            }
                            break;
                        case 2:
                            for (int x = 0; x <= edges_sub[e.id()]; x++) {
                                face_discr[f][x][edges_sub[ey]] = edge_discr[e.id()][x];
                                //std::cout << edge_discr[e.id()][x] << std::endl;
                            }
                            break;
                        default:
                            break;
                    }
                } else {
                    std::cout << "NOT ON FIRST EDGE : " << std::endl;
                    //std::cout << "\tEdge " << edge_discr[e.id()][0] << "," << edge_discr[e.id()][Ait]
                    //         << std::endl;
                }
            }
        }
        for (auto e: edges) {
            std::cout << "Edge " << edge_discr[e.id()][0] << "," << edge_discr[e.id()][edges_sub[e.id()]] << std::endl;
            if (nodes[3].id() == edge_discr[e.id()][0] && nodes[0].id() == edge_discr[e.id()][edges_sub[e.id()]]) {
                // + direction


                for (int y = 0; y <= edges_sub[e.id()]; y++) {
                    face_discr[f][0][y] = edge_discr[e.id()][edges_sub[e.id()] - y];
                    //std::cout << edge_discr[e.id()][edges_sub[e.id()] - y] << std::endl;
                }

            } else if (nodes[3].id() == edge_discr[e.id()][edges_sub[e.id()]] &&
                       nodes[0].id() == edge_discr[e.id()][0]) {
                // - direction
                for (int y = 0; y <= edges_sub[e.id()]; y++) {
                    face_discr[f][0][y] = edge_discr[e.id()][y];
                    //std::cout << edge_discr[e.id()][y] << std::endl;
                }

            } else {
                std::cout << "NOT ON FIRST EDGE : " << std::endl;
                std::cout << "\tEdge " << edge_discr[e.id()][0] << "," << edge_discr[e.id()][edges_sub[e.id()]]
                          << std::endl;
            }
        }
    }

    int edge_treated = m_hmesh->newMark<Edge>();

    std::set<TCellID> non_valid_edges;
    for(auto f : m_hmesh->faces()){
        if(linkerH_G.getGeomDim<Face>(f) == 3) {
            std::vector<cad::GeomSurface*> surfaces;
            manager.getSurfaces(surfaces);
            cad::GeomSurface* surface;
            for(auto s : surfaces){
                if(s->id() == linkerH_G.getGeomId<Face>(f)){
                    surface = s;
                }
            }
            std::vector<cad::GeomCurve*> curves;
            curves = surface->curves();


            Face face = m_hmesh->get<Face>(f);
            std::vector<TCellID> edges = face.getIDs<Edge>();
            for (auto e : edges) {
                if (linkerH_G.getGeomDim<Edge>(e) == 3 && !m_hmesh->isMarked<Edge>(e,edge_treated)){
                    m_hmesh->mark<Edge>(e,edge_treated);
                    std::cout<<"edge test "<<e<<std::endl;
                    bool non_valid_edge = false;
                    for(int i = 1; i<edges_sub[e]; i++) {
                        Node node = t_mesh.get<Node>(edge_discr[e][i]);
                        math::Point project_surf(node.getPoint().X(), node.getPoint().Y(), node.getPoint().Z());
                        surface->project(project_surf);
                        for(auto c : curves){
                            math::Point point(node.getPoint());
                            c->project(point);

                            //std::cout<<"distance "<<project_surf.distance(point)<<std::endl;
                            if (project_surf.distance(point) == 0) {
                                std::cout<<"curve "<<c->id()<<std::endl;
                                //std::cout<<"surf point "<<project_surf<<std::endl;
                                //std::cout<<"curv point "<<point<<std::endl;
                                non_valid_edge = true;
                            }
                        }
                    }
                    if(non_valid_edge) non_valid_edges.insert(e);
                }
            }
        }
    }

    for(auto e : non_valid_edges){
        if(linkerH_G.getGeomDim<Edge>(e) == 3) {
            Edge edge = m_hmesh->get<Edge>(e);

            math::Point p_first;
            math::Point p_second;

            if (edge.get<Node>()[0].id() < edge.get<Node>()[1].id()) {
                p_first = edge.get<Node>()[0].getPoint();
                p_second = edge.get<Node>()[1].getPoint();
            } else {
                p_first = edge.get<Node>()[1].getPoint();
                p_second = edge.get<Node>()[0].getPoint();
            }

            for (int i_p = 1; i_p < edges_sub[e]; i_p++) {
                //std::cout << "node on e " << edge_discr[e][i_p] << std::endl;
                math::Point new_p = ((1. - (1. / i_p)) * p_first) + ((1. / i_p) * p_second);
                std::vector<math::Point> points;
                points.push_back(t_mesh.get<Node>(edge_discr[e][i_p]).getPoint());
                for (auto f : edge.get<Face>()) {
                    if ((f.get<Region>().size() == 1 && linkerH_G.getGeomDim<Edge>(e) == 3) ||
                        (f.get<Region>().size() == 2 && linkerH_G.getGeomDim<Edge>(e) == 0)) {

                        int indice_i = 0;
                        int indice_j = 0;
                        for (int i = 0; i < face_discr[f.id()].size(); i++) {
                            for (int j = 0; j < face_discr[f.id()][i].size(); j++) {
                                //std::cout<<face_discr[f.id()][i][j]<<" "<<std::endl;
                                if (face_discr[f.id()][i][j] == edge_discr[e][i_p]) {
                                    //std::cout << "trouvé" << std::endl;
                                    indice_i = i;
                                    indice_j = j;
                                }
                            }
                            //std::cout<<std::endl;
                        }

                        int i_target = 0, j_target = 0;

                        if (indice_i == 0) {
                            i_target = face_discr[f.id()].size() - 1;
                            j_target = indice_j;
                        } else if (indice_i == face_discr[f.id()].size() - 1) {
                            i_target = 0;
                            j_target = indice_j;
                        } else if (indice_j == 0) {
                            i_target = indice_i;
                            j_target = face_discr[f.id()][0].size() - 1;
                        } else if (indice_j == face_discr[f.id()][0].size() - 1) {
                            i_target = indice_i;
                            j_target = 0;
                        }

                        points.push_back(t_mesh.get<Node>(face_discr[f.id()][i_target][j_target]).getPoint());


                    }
                }
                math::Point masscenter = math::Point::massCenter(points);
                //math::Point final_point = masscenter + new_p;
                //final_point.setXYZ(final_point.X(),final_point.Y(),final_point.Z());
                t_mesh.get<Node>(edge_discr[e][i_p]).setPoint(masscenter);
            }
        }
    }

    Variable<int>* face_id = t_mesh.getOrCreateVariable<int, GMDS_NODE>("face_id");

    for(auto f : m_hmesh->faces()){
        Face face = m_hmesh->get<Face>(f);
        int nb_div_i = face_discr[f].size()-1;
        int nb_div_j = face_discr[f][0].size()-1;
        //Transfinite interpolation to fill the interior of the face matrix
        math::Point p_first = face.get<Node>()[0].getPoint();
        math::Point p_second = face.get<Node>()[1].getPoint();
        math::Point p_ter = face.get<Node>()[2].getPoint();
        math::Point p_qua = face.get<Node>()[3].getPoint();


        for (double j = 1; j < nb_div_j; j++) {
            for (double i = 1; i < nb_div_i; i++) {

                math::Point point((p_second * (i / nb_div_i) + p_first * (1 - (i / nb_div_i))) * (1 - j / nb_div_j) +
                                  (p_ter * (i / nb_div_i) + p_qua * (1 - (i / nb_div_i))) * (j / nb_div_j));

                Node face_node = t_mesh.newNode(point);
                face_id->set(face_node.id(),f);
                //std::cout<<face_discr[f][int(i)][0]<<","<<face_discr[f][int(i)][int(iteration)]<<","<<face_discr[f][0][int(j)]<<","<<face_discr[f][int(iteration)][int(j)]<<std::endl;
                point = ((1-(i / nb_div_i))*(t_mesh.get<Node>(face_discr[f][0][int(j)]).getPoint()) + (i / nb_div_i)*(t_mesh.get<Node>(face_discr[f][int(nb_div_i)][int(j)]).getPoint())
                         + (1-(j / nb_div_j))*(t_mesh.get<Node>(face_discr[f][int(i)][0]).getPoint()) + (j / nb_div_j)*(t_mesh.get<Node>(face_discr[f][int(i)][int(nb_div_j)]).getPoint()))- point;
                face_discr[f][int(i)][int(j)] = face_node.id();
                face_node.setPoint(point);
                //std::cout<<"Point : "<<point<<std::endl;
                //m_hmesh->newTriangle(face_node,face_node,face_node);
                if(linkerH_G.getGeomDim<Face>(f) == 3){
                    linker_hexa.linkNodeToSurface(face_node.id(),linkerH_G.getGeomId<Face>(f));
                }
            }
        }

    }

    int nb_blocks = m_hmesh->getNbRegions();

    Variable<int>* block_id = t_mesh.newVariable<int, GMDS_REGION>("block_id");

    for(int r = 0; r < nb_blocks; r++){

        std::cout<<"start block"<<std::endl;

        Region block = m_hmesh->get<Region>(r);
        for(auto f : block.getIDs<Face>()){
            std::cout<<"face "<<f<<std::endl;
        }
        std::vector<Node> nodes = block.get<Node>();
        std::vector<Edge> edges = block.get<Edge>();
        std::vector<Face> faces = block.get<Face>();

        int subx;
        int suby;
        int subz;

        for(auto e : m_hmesh->edges()){
            Edge edge = m_hmesh->get<Edge>(e);
            if((edge.getIDs<Node>()[0] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[1] == block.getIDs<Node>()[1]) ||
                    (edge.getIDs<Node>()[1] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[0] == block.getIDs<Node>()[1])){
                subx = edges_sub[e];
            } else if((edge.getIDs<Node>()[0] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[1] == block.getIDs<Node>()[3]) ||
                      (edge.getIDs<Node>()[1] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[0] == block.getIDs<Node>()[3])){
                suby = edges_sub[e];
            } else if((edge.getIDs<Node>()[0] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[1] == block.getIDs<Node>()[4]) ||
                      (edge.getIDs<Node>()[1] == block.getIDs<Node>()[0] && edge.getIDs<Node>()[0] == block.getIDs<Node>()[4])){
                subz = edges_sub[e];
            }
        }

        TCellID block_discr[subx+1][suby+1][subz+1];


        Face face_0 = faces[0];
        int index_0 = 0;
        int index_1 = 0;
        std::vector<Node> f0_nodes = face_0.get<Node>();

        std::cout<<"Face block "<<nodes[0]<<","<<nodes[1]<<","<<nodes[2]<<","<<nodes[3]<<std::endl;
        std::cout<<"Face "<<f0_nodes[0]<<","<<f0_nodes[1]<<","<<f0_nodes[2]<<","<<f0_nodes[3]<<std::endl;

        for (int i = 0; i < 4; ++i) {
            std::cout<<f0_nodes[i].id()<<" ? "<<nodes[0].id()<<std::endl;
            if(f0_nodes[i].id()  == nodes[0].id()){
                index_0 = i;
                break;
            }
        }
        std::cout<<std::endl;
        for (int i = 0; i < 4; ++i) {
            std::cout<<f0_nodes[i].id()<<" ? "<<nodes[1].id()<<std::endl;
            if(f0_nodes[i].id() == nodes[1].id()){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;



        if(index_0 == 3 && index_1 == 0){
            for(int x = 0; x<subx+1; x++){
                for(int y = 0; y<suby+1; y++){
                    block_discr[x][y][0] = face_discr[face_0.id()][suby-y][x];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int x = 0; x<subx+1; x++){
                for(int y = 0; y<suby+1; y++){
                    block_discr[x][y][0] = face_discr[face_0.id()][y][x];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][subx-x][y];
                    }
                }
            }else if(index_1 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][suby-y][subx-x];
                    }
                }
            }else if(index_1 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][x][subx-y];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][x][y];
                    }
                }
            }else if(index_0 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][y][subx-x];
                    }
                }
            }else if(index_0 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][0] = face_discr[face_0.id()][subx-x][suby-y];
                    }
                }
            }
        }

        Face face_1 = faces[1];
        index_0 = 0;
        index_1 = 0;
        std::vector<Node> f1_nodes = face_1.get<Node>();

        std::cout<<"Face block "<<nodes[0].id()<<","<<nodes[1].id()<<","<<nodes[5].id()<<","<<nodes[4].id()<<std::endl;
        std::cout<<"Face "<<f1_nodes[0].id()<<","<<f1_nodes[1].id()<<","<<f1_nodes[2].id()<<","<<f1_nodes[3].id()<<std::endl;

        for (int i = 0; i < 4; ++i) {
            if(f1_nodes[i] == nodes[0]){
                index_0 = i;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if(f1_nodes[i] == nodes[1]){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;

        if(index_0 == 3 && index_1 == 0){
            for(int x = 0; x<subx+1; x++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[x][0][z] = face_discr[face_1.id()][subz-z][x];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int x = 0; x<subx+1; x++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[x][0][z] = face_discr[face_1.id()][z][x];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][subx-x][z];
                    }
                }
            }else if(index_1 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][subz-z][subx-x];
                    }
                }
            }else if(index_1 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][x][subz-z];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][x][z];
                    }
                }
            }else if(index_0 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][z][subx-x];
                    }
                }
            }else if(index_0 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][0][z] = face_discr[face_1.id()][subx-x][subz-z];
                    }
                }
            }
        }

        Face face_2 = faces[2];
        index_0 = 0;
        index_1 = 0;
        std::vector<Node> f2_nodes = face_2.get<Node>();

        std::cout<<"Face block "<<nodes[1].id()<<","<<nodes[2].id()<<","<<nodes[6].id()<<","<<nodes[5].id()<<std::endl;
        std::cout<<"Face "<<f2_nodes[0].id()<<","<<f2_nodes[1].id()<<","<<f2_nodes[2].id()<<","<<f2_nodes[3].id()<<std::endl;

        for (int i = 0; i < 4; ++i) {
            if(f2_nodes[i] == nodes[1]){
                index_0 = i;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if(f2_nodes[i] == nodes[2]){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;


        if(index_0 == 3 && index_1 == 0){
            for(int y = 0; y<suby+1; y++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[subx][y][z] = face_discr[face_2.id()][subz-z][y];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int y = 0; y<suby+1; y++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[subx][y][z] = face_discr[face_2.id()][z][y];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][suby-y][z];
                    }
                }
            }else if(index_1 == 1){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][subz-z][suby-y];
                    }
                }
            }else if(index_1 == 2){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][y][subz-z];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][y][z];
                    }
                }
            }else if(index_0 == 1){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][z][suby-y];
                    }
                }
            }else if(index_0 == 2){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[subx][y][z] = face_discr[face_2.id()][suby-y][subz-z];
                    }
                }
            }
        }

        Face face_3 = faces[3];
        index_0 = 0;
        index_1 = 0;
        std::vector<Node> f3_nodes = face_3.get<Node>();

        std::cout<<"Face block "<<nodes[2].id()<<","<<nodes[3].id()<<","<<nodes[7].id()<<","<<nodes[6].id()<<std::endl;
        std::cout<<"Face "<<f3_nodes[0].id()<<","<<f3_nodes[1].id()<<","<<f3_nodes[2].id()<<","<<f3_nodes[3].id()<<std::endl;

        for (int i = 0; i < 4; ++i) {
            if(f3_nodes[i] == nodes[2]){
                index_0 = i;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if(f3_nodes[i] == nodes[3]){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;


        if(index_0 == 3 && index_1 == 0){
            for(int x = 0; x<subx+1; x++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[x][suby][z] = face_discr[face_3.id()][z][x];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int x = 0; x<subx+1; x++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[x][suby][z] = face_discr[face_3.id()][subz-z][x];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][x][z];
                    }
                }
            }else if(index_1 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][z][subx-x];
                    }
                }
            }else if(index_1 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][subx-x][subz-z];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][subx-x][z];
                    }
                }
            }else if(index_0 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][subz-z][subx-x];
                    }
                }
            }else if(index_0 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[x][suby][z] = face_discr[face_3.id()][x][subz-z];
                    }
                }
            }
        }

        Face face_4 = faces[4];
        index_0 = 0;
        index_1 = 0;
        std::vector<Node> f4_nodes = face_4.get<Node>();

        std::cout<<"Face block "<<nodes[3].id()<<","<<nodes[0].id()<<","<<nodes[4].id()<<","<<nodes[7].id()<<std::endl;
        std::cout<<"Face "<<f4_nodes[0].id()<<","<<f4_nodes[1].id()<<","<<f4_nodes[2].id()<<","<<f4_nodes[3].id()<<std::endl;

        for (int i = 0; i < 4; ++i) {
            if(f4_nodes[i] == nodes[3]){
                index_0 = i;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if(f4_nodes[i] == nodes[0]){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;

        std::cout<<block_discr[0][0][0]<<","<<block_discr[subx][0][0]<<","<<block_discr[subx][suby][0]<<","<<block_discr[0][suby][0]<<","<<
                 block_discr[0][0][subz]<<","<<block_discr[subx][0][subz]<<","<<block_discr[subx][suby][subz]<<","<<block_discr[0][suby][subz]<<std::endl;


        if(index_0 == 3 && index_1 == 0){
            for(int y = 0; y<suby+1; y++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[0][y][z] = face_discr[face_4.id()][z][y];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int y = 0; y<suby+1; y++){
                for(int z = 0; z<subz+1; z++){
                    block_discr[0][y][z] = face_discr[face_4.id()][subz-z][y];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][y][z];
                    }
                }
            }else if(index_1 == 1){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][z][suby-y];
                    }
                }
            }else if(index_1 == 2){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][suby-y][subz-z];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][suby-y][z];
                    }
                }
            }else if(index_0 == 1){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][subz-z][suby-y];
                    }
                }
            }else if(index_0 == 2){
                for(int y = 0; y<suby+1; y++){
                    for(int z = 0; z<subz+1; z++){
                        block_discr[0][y][z] = face_discr[face_4.id()][y][subz-z];
                    }
                }
            }
        }
        std::cout<<"after face 4"<<std::endl;

        Face face_5 = faces[5];
        index_0 = 0;
        index_1 = 0;
        std::vector<Node> f5_nodes = face_5.get<Node>();

        std::cout<<"Face block "<<nodes[4]<<","<<nodes[5]<<","<<nodes[6]<<","<<nodes[7]<<std::endl;
        std::cout<<"Face "<<f5_nodes[0]<<","<<f5_nodes[1]<<","<<f5_nodes[2]<<","<<f5_nodes[3]<<std::endl;

        for (int i = 0; i < 4; ++i) {
            if(f5_nodes[i] == nodes[4]){
                index_0 = i;
                break;
            }
        }
        for (int i = 0; i < 4; ++i) {
            if(f5_nodes[i] == nodes[5]){
                index_1 = i;
                break;
            }
        }

        std::cout<<"index_0 = "<<index_0<<std::endl;
        std::cout<<"index_1 = "<<index_1<<std::endl;


        if(index_0 == 3 && index_1 == 0){
            for(int x = 0; x<subx+1; x++){
                for(int y = 0; y<suby+1; y++){
                    block_discr[x][y][subz] = face_discr[face_5.id()][suby-y][x];
                }
            }
        }else if(index_0 == 0 && index_1 == 3){
            for(int x = 0; x<subx+1; x++){
                for(int y = 0; y<suby+1; y++){
                    block_discr[x][y][subz] = face_discr[face_5.id()][y][x];
                }
            }
        }else if(index_0 == index_1+1){
            if(index_1 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][subx-x][y];
                    }
                }
            }else if(index_1 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][suby-y][subx-x];
                    }
                }
            }else if(index_1 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][x][suby-y];
                    }
                }
            }
        }else{
            if(index_0 == 0){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][x][y];
                    }
                }
            }else if(index_0 == 1){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][y][subx-x];
                    }
                }
            }else if(index_0 == 2){
                for(int x = 0; x<subx+1; x++){
                    for(int y = 0; y<suby+1; y++){
                        block_discr[x][y][subz] = face_discr[face_5.id()][subx-x][suby-y];
                    }
                }
            }
        }

        for (double z = 1; z <= subz-1; z++) {
            for (double y = 1; y <= suby-1; y++) {
                for (double x = 1; x <= subx-1; x++) {

                    math::Point p_new;

                    math::Point pf0 = t_mesh.get<Node>(block_discr[int(x)][int(y)][0]).getPoint();
                    math::Point pf1 = t_mesh.get<Node>(block_discr[int(x)][0][int(z)]).getPoint();
                    math::Point pf2 = t_mesh.get<Node>(block_discr[subx][int(y)][int(z)]).getPoint();
                    math::Point pf3 = t_mesh.get<Node>(block_discr[int(x)][suby][int(z)]).getPoint();
                    math::Point pf4 = t_mesh.get<Node>(block_discr[0][int(y)][int(z)]).getPoint();
                    math::Point pf5 = t_mesh.get<Node>(block_discr[int(x)][int(y)][subz]).getPoint();

                    math::Point p0 = nodes[0].getPoint();
                    math::Point p1 = nodes[1].getPoint();
                    math::Point p2 = nodes[2].getPoint();
                    math::Point p3 = nodes[3].getPoint();
                    math::Point p4 = nodes[4].getPoint();
                    math::Point p5 = nodes[5].getPoint();
                    math::Point p6 = nodes[6].getPoint();
                    math::Point p7 = nodes[7].getPoint();

                    std::vector<math::Point> points;
                    points.push_back(pf0);
                    points.push_back(pf1);
                    points.push_back(pf2);
                    points.push_back(pf3);
                    points.push_back(pf4);
                    points.push_back(pf5);

                    p_new = math::Point::massCenter(points);


                    /*p_new = /*(1 - (z / iteration))*pf0 + (z / iteration)*pf5 +
                            (1 - (y / iteration))*pf1 + (y / iteration)*pf3 +
                            (1 - (x / iteration))*pf4 + (x / iteration)*pf2 -
                            (((p0 * (1 - x / iteration) + p1 * ((x / iteration))) * (1 -y / iteration) +
                              (p3 * (1 - x / iteration) + p2 * ((x / iteration))) * ( y / iteration)) * (1 -z / iteration) +
                             ((p4 * (1 - x / iteration) + p5 * ((x / iteration))) * (1 -y / iteration) +
                              (p7 * (1 - x / iteration) + p6 * ((x / iteration))) * ( y / iteration)) * (z / iteration));
*/
                    //std::cout<<p_new<<std::endl;

                    Node node_new = t_mesh.newNode(p_new);
                    block_discr[int(x)][int(y)][int(z)] = node_new.id();
                    //std::cout<<"test"<<std::endl;
                    //m_hmesh->newTriangle(node_new.id(), node_new.id(), node_new.id());
                }
            }
        }

        for (int z = 0; z < subz; z++) {
            for (int y = 0; y < suby; y++) {
                for (int x = 0; x < subx; x++) {
                    Region hex = t_mesh.newHex(block_discr[x][y][z],block_discr[x+1][y][z],block_discr[x+1][y+1][z],block_discr[x][y+1][z],
                                               block_discr[x][y][z+1],block_discr[x+1][y][z+1],block_discr[x+1][y+1][z+1],block_discr[x][y+1][z+1]);
                    //(*m_block_id)[hex.id()] = r;

                    //std::cout<<block_discr[x][y][z]<<","<<block_discr[x+1][y][z]<<","<<block_discr[x+1][y+1][z]<<","<<block_discr[x][y+1][z]<<","<<
                    //         block_discr[x][y][z+1]<<","<<block_discr[x+1][y][z+1]<<","<<block_discr[x+1][y+1][z+1]<<","<<block_discr[x][y+1][z+1]<<std::endl;

                    block_id->set(hex.id(), r);
                }
            }
        }
        //break;

    }



    std::cout<<"Nb blocks "<<blocks.size()<<std::endl;

    for(auto b:blocks){
        std::vector<Node> b_nodes = b.get<Node>();
        for(auto n:b_nodes){
            //n.remove(b);
        }
        //m_hmesh->deleteRegion(b);
    }

    tdoc.buildFacesAndR2F();
    tdoc.buildEdgesAndX2E();
    tdoc.updateUpwardConnectivity();


    for(auto n : t_mesh.nodes()){
        Node node = t_mesh.get<Node>(n);
        if(linker_hexa.getGeomDim<Node>(n) == 2){
            for(auto e : node.get<Edge>()){
                if(e.getIDs<Node>()[0] == n && (linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[1]) > 0 && linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[1]) < 3)){
                    linker_hexa.linkToCurve(e,linker_hexa.getGeomId<Node>(n));
                }else if(e.getIDs<Node>()[1] == n && (linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[0]) > 0 && linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[0]) < 3)){
                    linker_hexa.linkToCurve(e,linker_hexa.getGeomId<Node>(n));
                }
            }
        }
        else if(linker_hexa.getGeomDim<Node>(n) == 3){
            for(auto e : node.get<Edge>()){
                if(e.getIDs<Node>()[0] == n && (linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[1]) > 0)){
                    linker_hexa.linkToSurface(e,linker_hexa.getGeomId<Node>(n));
                }else if(e.getIDs<Node>()[1] == n && (linker_hexa.getGeomDim<Node>(e.getIDs<Node>()[0]) > 0)){
                    linker_hexa.linkToSurface(e,linker_hexa.getGeomId<Node>(n));
                }
            }
        }
    }
    for(auto f : t_mesh.faces()){
        Face face = t_mesh.get<Face>(f);
        bool on_surf = true;
        int geom_id = 0;
        for(auto n : face.getIDs<Node>()){
            if(linker_hexa.getGeomDim<Node>(n) == 0){
                on_surf = false;
            }
            if(linker_hexa.getGeomDim<Node>(n) == 3){
                geom_id = linker_hexa.getGeomId<Node>(n);
            }
        }
        if(on_surf) linker_hexa.linkFaceToSurface(f,geom_id);
    }


    Variable<int>* n_on_point = t_mesh.getOrCreateVariable<int,GMDS_NODE>("n_on_point");
    Variable<int>* n_on_curve = t_mesh.getOrCreateVariable<int,GMDS_NODE>("n_on_curve");
    Variable<int>* n_on_surface = t_mesh.getOrCreateVariable<int,GMDS_NODE>("n_on_surface");
    for(auto n : t_mesh.nodes()){
        n_on_point->set(n,-1);
        n_on_curve->set(n,-1);
        n_on_surface->set(n,-1);
        if(linker_hexa.getGeomDim<Node>(n) == 1){
            n_on_point->set(n, linker_hexa.getGeomId<Node>(n));
        }else if(linker_hexa.getGeomDim<Node>(n) == 2){
            n_on_curve->set(n, linker_hexa.getGeomId<Node>(n));
        }else if(linker_hexa.getGeomDim<Node>(n) == 3){
            n_on_surface->set(n, linker_hexa.getGeomId<Node>(n));
        }
    }


    IGMeshIOService link_ioService(linker_hexa.mesh());
    VTKWriter linkvtkWriter(&link_ioService);
    linkvtkWriter.setCellOptions(gmds::N|gmds::R);
    linkvtkWriter.setDataOptions(gmds::N|gmds::R);
    linkvtkWriter.write("/home/simon/Data/Results_debug/linker.vtk");

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write("/home/simon/Data/Results_debug/inter.vtk");

    cad::GeomSmoother smoother(&linker_hexa);
    smoother.smoothCurves(10);
    std::cout<<"curves"<<std::endl;
    smoother.smoothSurfaces(10);
    std::cout<<"surfaces"<<std::endl;
    smoother.smoothVolumes(5);
    std::cout<<"volumes"<<std::endl;

    IGMeshIOService t_ioService(&t_mesh);
    VTKWriter tvtkWriter(&t_ioService);
    tvtkWriter.setCellOptions(gmds::N|gmds::R);
    tvtkWriter.setDataOptions(gmds::N|gmds::R);
    tvtkWriter.write("/home/simon/Data/Results_debug/t_inter.vtk");

}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::interpolation(int Ait) {

    std::cout<<"Interpolation"<<std::endl;

    double iteration = Ait;

    for(auto r : m_hmesh->regions()) {

        TCellID tab[int(iteration)+1][int(iteration)+1][int(iteration)+1];

        std::cout<<"Block : "<<r<<std::endl;
        std::vector<Node> nodes = m_hmesh->get<Region>(r).get<Node>();
        for (auto n : nodes) {
            std::cout << n.id() << ",";
        }
        //std::cout << std::endl;
        std::vector<Edge> edges_ = m_hmesh->get<Region>(r).get<Edge>();
        for (auto e : edges_) {
            //std::cout << e.get<Node>()[0]<<", "<< e.get<Node>()[1]<<std::endl;
        }
        std::cout << std::endl;

        tab[0][0][0] = nodes[0].id();
        math::Point p0 = nodes[0].getPoint();

        //DEBUG
        (*m_X)[nodes[0].id()] = 0;
        (*m_Y)[nodes[0].id()] = 0;
        (*m_Z)[nodes[0].id()] = 0;

        tab[int(iteration)][0][0] = nodes[1].id();
        math::Point p1 = nodes[1].getPoint();

        //DEBUG
        (*m_X)[nodes[1].id()] = int(iteration);
        (*m_Y)[nodes[1].id()] = 0;
        (*m_Z)[nodes[1].id()] = 0;

        tab[int(iteration)][int(iteration)][0] = nodes[2].id();
        math::Point p2 = nodes[2].getPoint();

        //DEBUG
        (*m_X)[nodes[2].id()] = int(iteration);
        (*m_Y)[nodes[2].id()] = int(iteration);
        (*m_Z)[nodes[2].id()] = 0;

        tab[0][int(iteration)][0] = nodes[3].id();
        math::Point p3 = nodes[3].getPoint();

        //DEBUG
        (*m_X)[nodes[3].id()] = 0;
        (*m_Y)[nodes[3].id()] = int(iteration);
        (*m_Z)[nodes[3].id()] = 0;

        tab[0][0][int(iteration)] = nodes[4].id();
        math::Point p4 = nodes[4].getPoint();

        //DEBUG
        (*m_X)[nodes[4].id()] = 0;
        (*m_Y)[nodes[4].id()] = 0;
        (*m_Z)[nodes[4].id()] = int(iteration);

        tab[int(iteration)][0][int(iteration)] = nodes[5].id();
        math::Point p5 = nodes[5].getPoint();

        //DEBUG
        (*m_X)[nodes[5].id()] = int(iteration);
        (*m_Y)[nodes[5].id()] = 0;
        (*m_Z)[nodes[5].id()] = int(iteration);

        tab[int(iteration)][int(iteration)][int(iteration)] = nodes[6].id();
        math::Point p6 = nodes[6].getPoint();

        //DEBUG
        (*m_X)[nodes[6].id()] = int(iteration);
        (*m_Y)[nodes[6].id()] = int(iteration);
        (*m_Z)[nodes[6].id()] = int(iteration);

        tab[0][int(iteration)][int(iteration)] = nodes[7].id();
        math::Point p7 = nodes[7].getPoint();

        //DEBUG
        (*m_X)[nodes[7].id()] = 0;
        (*m_Y)[nodes[7].id()] = int(iteration);
        (*m_Z)[nodes[7].id()] = int(iteration);


        math::Point p_new;

        std::vector<std::pair<Node, Node>> edges;

        edges.emplace_back(nodes[0], nodes[1]);
        edges.emplace_back(nodes[1], nodes[2]);
        edges.emplace_back(nodes[3], nodes[2]);
        edges.emplace_back(nodes[0], nodes[3]);
        edges.emplace_back(nodes[4], nodes[5]);
        edges.emplace_back(nodes[5], nodes[6]);
        edges.emplace_back(nodes[7], nodes[6]);
        edges.emplace_back(nodes[4], nodes[7]);
        edges.emplace_back(nodes[0], nodes[4]);
        edges.emplace_back(nodes[1], nodes[5]);
        edges.emplace_back(nodes[2], nodes[6]);
        edges.emplace_back(nodes[3], nodes[7]);


        //--------------------------------------------------
        //Edges
        int it_edge = 0;
        int type_edges[12];

        for (auto edge : edges) {

            int geom_dim_e1,geom_dim_e2;
            int geom_id_e1,geom_id_e2;

            std::vector<Edge> r_edges = m_hmesh->get<Region>(r).get<Edge>();
            Edge current_r_edge;
            bool border_edge = false;
            if((linkerH_G.getGeomDim(edge.first) == 2 || linkerH_G.getGeomDim(edge.first) == 3) &&
               ((linkerH_G.getGeomDim(edge.second) == 2 || linkerH_G.getGeomDim(edge.second) == 3))) {
                for (auto e : r_edges) {
                    if ((e.get<Node>()[0] == edge.first && e.get<Node>()[1] == edge.second) ||
                        (e.get<Node>()[1] == edge.first && e.get<Node>()[0] == edge.second)) {
                        current_r_edge = e;
                        for (auto f : e.get<Face>()) {
                            //std::cout<<"test"<<std::endl;
                            if (f.get<Region>().size() == 1) {
                                border_edge = true;
                            }
                        }
                    }
                }
            }else{
                for (auto e : r_edges) {
                    if ((e.get<Node>()[0] == edge.first && e.get<Node>()[1] == edge.second) ||
                        (e.get<Node>()[1] == edge.first && e.get<Node>()[0] == edge.second)) {
                        current_r_edge = e;
                    }
                }
            }

            int geom_id = -1;
            int type = 0;

            if(border_edge){

                geom_id_e1 = linkerH_G.getGeomId(edge.first);
                geom_id_e2 = linkerH_G.getGeomId(edge.second);

                if(linkerH_G.getGeomDim(edge.first) == 2 && linkerH_G.getGeomDim(edge.second) == 2){
                    if(geom_id_e1 != geom_id_e2) {
                        int gc1 = geom_id_e1;
                        int gc2 = geom_id_e2;

                        std::vector<cad::GeomCurve *> curves;
                        manager.getCurves(curves);

                        cad::GeomCurve *GC1;
                        cad::GeomCurve *GC2;

                        for (auto c : curves) {
                            if (c->id() == gc1) {
                                GC1 = c;
                            }
                            if (c->id() == gc2) {
                                GC2 = c;
                            }
                        }
                        geom_id = manager.getCommonSurface(GC1, GC2);
                        type = 3;
                        linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                    }else{
                        geom_id = geom_id_e1;
                        type = 2;
                        linkerH_G.linkEdgeToCurve(current_r_edge.id(),geom_id);
                    }
                }else if(linkerH_G.getGeomDim(edge.first) == 3 && linkerH_G.getGeomDim(edge.second) == 3){
                    if(geom_id_e1 != geom_id_e2) {
                        type = 4;
                    } else{
                        geom_id = geom_id_e1;
                        type = 3;
                        linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                    }
                }else if(linkerH_G.getGeomDim(edge.first) == 3){

                    geom_id = geom_id_e1;
                    type = 3;
                    linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                }else if(linkerH_G.getGeomDim(edge.second) == 3){

                    geom_id = geom_id_e2;
                    type = 3;
                    linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                }

            }
            else {
                if (linkerH_G.getGeomDim(edge.first) <= linkerH_G.getGeomDim(edge.second)) {
                    geom_dim_e1 = linkerH_G.getGeomDim(edge.first);
                    geom_dim_e2 = linkerH_G.getGeomDim(edge.second);

                    geom_id_e1 = linkerH_G.getGeomId(edge.first);
                    geom_id_e2 = linkerH_G.getGeomId(edge.second);
                } else {
                    geom_dim_e1 = linkerH_G.getGeomDim(edge.second);
                    geom_dim_e2 = linkerH_G.getGeomDim(edge.first);

                    geom_id_e1 = linkerH_G.getGeomId(edge.second);
                    geom_id_e2 = linkerH_G.getGeomId(edge.first);
                }

                if (geom_dim_e1 == 1) {
                    if (geom_dim_e2 == 1) {
                        int gp1 = geom_id_e1;
                        int gp2 = geom_id_e2;

                        std::vector<cad::GeomPoint *> points;
                        manager.getPoints(points);

                        cad::GeomPoint *GP1;
                        cad::GeomPoint *GP2;

                        for (auto p : points) {
                            if (p->id() == gp1) {
                                GP1 = p;
                            }
                            if (p->id() == gp2) {
                                GP2 = p;
                            }
                        }
                        geom_id = manager.getCommonCurve(GP1, GP2);
                        type = 2;
                        linkerH_G.linkEdgeToCurve(current_r_edge.id(),geom_id);
                    } else if (geom_dim_e2 == 2) {

                        int gp1 = geom_id_e1;
                        std::vector<cad::GeomPoint *> points;
                        manager.getPoints(points);

                        cad::GeomPoint *GP1;

                        for (auto p : points) {
                            if (p->id() == gp1) {
                                GP1 = p;
                            }
                        }

                        math::Point point_gp1 = GP1->point();

                        std::vector<cad::GeomCurve *> curves;
                        manager.getCurves(curves);

                        cad::GeomCurve *GC1;

                        for (auto c : curves) {
                            if (c->id() == geom_id_e2) {
                                GC1 = c;
                            }
                        }

                        math::Point closest = GC1->closestPoint(point_gp1);

                        std::cout << "point_gp1 = " << point_gp1 << std::endl;
                        std::cout << "closest = " << closest << std::endl;

                        if (closest.X() != point_gp1.X() || closest.Y() != point_gp1.Y() ||
                            closest.Z() != point_gp1.Z()) {
                            std::cout << "point---edge" << std::endl;
                            type = 4;
                        } else {
                            geom_id = geom_id_e2;
                            type = 2;
                            linkerH_G.linkEdgeToCurve(current_r_edge.id(),geom_id);
                        }

                    } else if (geom_dim_e2 == 3) {
                        geom_id = geom_id_e2;
                        type = 3;
                        linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                    } else {
                        type = 4;
                    }
                } else if (geom_dim_e1 == 2) {
                    if (geom_dim_e2 == 2) {
                        if (geom_id_e1 == geom_id_e2) {
                            geom_id = geom_id_e1;
                            type = 2;
                            linkerH_G.linkEdgeToCurve(current_r_edge.id(),geom_id);
                        } else {
                            int gc1 = geom_id_e1;
                            int gc2 = geom_id_e2;

                            std::vector<cad::GeomCurve *> curves;
                            manager.getCurves(curves);

                            cad::GeomCurve *GC1;
                            cad::GeomCurve *GC2;

                            for (auto c : curves) {
                                if (c->id() == gc1) {
                                    GC1 = c;
                                }
                                if (c->id() == gc2) {
                                    GC2 = c;
                                }
                            }

                            geom_id = manager.getCommonSurface(GC1, GC2);

                            if(geom_id==-1){
                                type = 4;
                            }
                            else{
                                type = 3;
                                linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                            }
                        }
                    } else if (geom_dim_e2 == 3) {
                        type = 4;
                    } else {
                        type = 4;
                    }
                } else if (geom_dim_e1 == 3) {
                    if (geom_dim_e2 == 3) {
                        if (geom_id_e1 == geom_id_e2) {
                            geom_id = geom_id_e1;
                            type = 3;
                            linkerH_G.linkEdgeToSurface(current_r_edge.id(),geom_id);
                        }
                    } else {
                        type = 4;
                    }
                }
            }
            type_edges[it_edge] = type;

            math::Point p_first = edge.first.getPoint();
            math::Point p_second = edge.second.getPoint();

            math::Vector3d vector(p_first, p_second);
            math::Vector3d mini_vec = vector / 10;

            for (double i = 1; i < iteration; i++) {
                //math::Point point(p1+mini_vec*i);
                math::Point point;
                point = (p_second * (i / iteration) + p_first * (1 - (i / iteration)));

                //std::cout<<"point:"<<point<<std::endl;
                //std::cout<<"vector Y: "<<vector.Y()<<std::endl;
                if (type == 2) {
                    std::vector<cad::GeomCurve *> curves;
                    manager.getCurves(curves);
                    for (auto c : curves) {
                        if (c->id() == geom_id) {
                            c->project(point);
                        }
                    }
                } else if (type == 3) {
                    std::vector<cad::GeomSurface *> surfs;
                    manager.getSurfaces(surfs);
                    for (auto s : surfs) {
                        if (s->id() == geom_id) {
                            s->project(point);
                        }
                    }
                }
                Node inter_node = m_hmesh->newNode(point);
                if (type == 2) {
                    linkerH_G.linkToCurve(inter_node,geom_id);
                } else if (type == 3) {
                    linkerH_G.linkToSurface(inter_node,geom_id);
                }

                switch (it_edge) {
                    case 0 :
                        tab[int(i)][0][0] = inter_node.id();
                        break;
                    case 1 :
                        tab[int(iteration)][int(i)][0] = inter_node.id();
                        break;
                    case 2 :
                        tab[int(i)][int(iteration)][0] = inter_node.id();
                        break;
                    case 3 :
                        tab[0][int(i)][0] = inter_node.id();
                        break;
                    case 4 :
                        tab[int(i)][0][int(iteration)] = inter_node.id();
                        break;
                    case 5 :
                        tab[int(iteration)][int(i)][int(iteration)] = inter_node.id();
                        break;
                    case 6 :
                        tab[int(i)][int(iteration)][int(iteration)] = inter_node.id();
                        break;
                    case 7 :
                        tab[0][int(i)][int(iteration)] = inter_node.id();
                        break;
                    case 8 :
                        tab[0][0][int(i)] = inter_node.id();
                        break;
                    case 9 :
                        tab[int(iteration)][0][int(i)] = inter_node.id();
                        break;
                    case 10 :
                        tab[int(iteration)][int(iteration)][int(i)] = inter_node.id();
                        break;
                    case 11 :
                        tab[0][int(iteration)][int(i)] = inter_node.id();
                        break;
                    default:
                        break;
                }

                //m_hmesh->newTriangle(inter_node.id(), inter_node.id(), inter_node.id());
            }
            it_edge++;
        }



        //--------------------------------------------------
        //Faces

        int cpt = 0;

        int it_faces = 0;

        std::vector<Face> faces = m_hmesh->get<Region>(r).get<Face>();

        for (auto face : faces) {

            //std::cout<<"it_faces = "<<it_faces<<std::endl;
            int geom_id = -1;
            int type = 0;

            if(face.get<Region>().size() == 1){
                type = 3;
                for(auto e : face.get<Edge>()){
                    if(linkerH_G.getGeomDim(e) == 3){
                        geom_id = linkerH_G.getGeomId(e);
                    }
                }
            }else{
                type = 4;
            }
            std::cout<<"FACE : "<<face.getIDs<Node>()[0]<<","<<face.getIDs<Node>()[1]<<","<<face.getIDs<Node>()[2]<<","<<face.getIDs<Node>()[3]<<" Type : "<<type<<std::endl;


            math::Point p_first = face.get<Node>()[0].getPoint();
            math::Point p_second = face.get<Node>()[1].getPoint();
            math::Point p_ter = face.get<Node>()[2].getPoint();
            math::Point p_qua = face.get<Node>()[3].getPoint();


            for (double j = 1; j < iteration; j++) {
                for (double i = 1; i < iteration; i++) {
                    //math::Point point(p1+mini_vec*i);
                    math::Point point((p_second * (i / iteration) + p_first * (1 - (i / iteration))) * (1 - j / iteration) +
                                      (p_ter * (i / iteration) + p_qua * (1 - (i / iteration))) * (j / iteration));

                    //std::cout<<"point:"<<point<<std::endl;
                    //std::cout<<"vector Y: "<<vector.Y()<<std::endl;

                    Node face_node = m_hmesh->newNode(point);
                    if(type == 3){
                        linkerH_G.linkToSurface(face_node,geom_id);
                    }
                    switch (it_faces) {
                        case 0 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[int(j)][0][0]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[int(j)][int(iteration)][0]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[0][int(i)][0]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[int(iteration)][int(i)][0]).getPoint()))- point;
                            tab[int(j)][int(i)][0] = face_node.id();
                            face_node.setPoint(point);

                            break;
                        case 1 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[0][0][int(j)]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[int(iteration)][0][int(j)]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[int(i)][0][0]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[int(i)][0][int(iteration)]).getPoint()))- point;
                            tab[int(i)][0][int(j)] = face_node.id();
                            face_node.setPoint(point);
                            break;
                        case 2 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[int(iteration)][0][int(j)]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[int(iteration)][int(iteration)][int(j)]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[int(iteration)][int(i)][0]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[int(iteration)][int(i)][int(iteration)]).getPoint()))- point;
                            tab[int(iteration)][int(i)][int(j)] = face_node.id();
                            face_node.setPoint(point);
                            break;
                        case 3 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[int(iteration)][int(iteration)][int(j)]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[0][int(iteration)][int(j)]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[int(iteration)-int(i)][int(iteration)][0]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[int(iteration)-int(i)][int(iteration)][int(iteration)]).getPoint()))- point;
                            tab[int(iteration)-int(i)][int(iteration)][int(j)] = face_node.id();
                            face_node.setPoint(point);
                            break;
                        case 4 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[0][int(j)][0]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[0][int(j)][int(iteration)]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[0][0][int(i)]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[0][int(iteration)][int(i)]).getPoint()))- point;
                            tab[0][int(j)][int(i)] = face_node.id();
                            face_node.setPoint(point);
                            break;
                        case 5 :
                            point = ((1-(i / iteration))*(m_hmesh->get<Node>(tab[0][int(j)][int(iteration)]).getPoint()) + (i / iteration)*(m_hmesh->get<Node>(tab[int(iteration)][int(j)][int(iteration)]).getPoint())
                                     + (1-(j / iteration))*(m_hmesh->get<Node>(tab[int(i)][0][int(iteration)]).getPoint()) + (j / iteration)*(m_hmesh->get<Node>(tab[int(i)][int(iteration)][int(iteration)]).getPoint()))- point;
                            tab[int(i)][int(j)][int(iteration)] = face_node.id();
                            face_node.setPoint(point);
                            break;

                        default:
                            break;
                            //std::cout<<"default face "<<std::endl;
                    }
                    if (type == 3) {
                        std::vector<cad::GeomSurface *> surfs;
                        manager.getSurfaces(surfs);
                        for (auto s : surfs) {
                            if (s->id() == geom_id) {
                                //s->project(point);
                            }
                        }
                    }
                    //face_node.setPoint(point);


                    //m_hmesh->newTriangle(face_node.id(), face_node.id(), face_node.id());
                }
            }
            it_faces++;
        }

        //--------------------------------------------------
        //Volume

        for (double z = 1; z <= iteration-1; z++) {
            for (double y = 1; y <= iteration-1; y++) {
                for (double x = 1; x <= iteration-1; x++) {
                    p_new = ((p0 * (1 - x / iteration) + p1 * ((x / iteration))) * (1 -y / iteration) +
                             (p3 * (1 - x / iteration) + p2 * ((x / iteration))) * ( y / iteration)) * (1 -z / iteration) +
                            ((p4 * (1 - x / iteration) + p5 * ((x / iteration))) * (1 -y / iteration) +
                             (p7 * (1 - x / iteration) + p6 * ((x / iteration))) * ( y / iteration)) * (z / iteration);

                    Node node_new = m_hmesh->newNode(p_new);
                    tab[int(x)][int(y)][int(z)] = node_new.id();
                    //std::cout<<"test"<<std::endl;
                    //m_hmesh->newTriangle(node_new.id(), node_new.id(), node_new.id());
                }
            }
        }
        for (int z = 0; z < iteration+1; z++) {
            for (int y = 0; y < iteration+1; y++) {
                for (int x = 0; x < iteration+1; x++) {

                    (*m_X)[tab[x][y][z]] = x;
                    (*m_Y)[tab[x][y][z]] = y;
                    (*m_Z)[tab[x][y][z]] = z;
                }
            }
        }

        for (int z = 0; z < iteration; z++) {
            for (int y = 0; y < iteration; y++) {
                for (int x = 0; x < iteration; x++) {
                    Region hex = m_hmesh->newHex(tab[x][y][z],tab[x+1][y][z],tab[x+1][y+1][z],tab[x][y+1][z],
                                                 tab[x][y][z+1],tab[x+1][y][z+1],tab[x+1][y+1][z+1],tab[x][y+1][z+1]);
                    (*m_block_id)[hex.id()] = r;
                }
            }
        }


        IGMeshIOService ioService(m_hmesh);
        VTKWriter vtkWriter(&ioService);
        vtkWriter.setCellOptions(gmds::N|gmds::R);
        vtkWriter.setDataOptions(gmds::N|gmds::R);
        vtkWriter.write("/ccc/temp/cont001/ocre/calderans/Results_debug/inter.vtk");

        if(r == 4)
            //break;

            for (int z = 0; z < iteration+1; z++) {
                for (int y = 0; y < iteration+1; y++) {
                    for (int x = 0; x < iteration+1; x++) {

                        (*m_X)[tab[x][y][z]] = 0;
                        (*m_Y)[tab[x][y][z]] = 0;
                        (*m_Z)[tab[x][y][z]] = 0;
                    }
                }
            }
    }


}
/*----------------------------------------------------------------------------*/

void DualBlockingSession::refineSurfaceSheet() {


    for(int i = 0; i<dual_sheets.size();i++){

        if(dual_sheets[i].isBoundary() /*&& dual_sheets[i].getID() == ASheetID*/) {

            int id_sheet = dual_sheets[i].getID();
            std::vector<TCellID> sheet_tets = dual_sheets[i].getSurface();

            for (auto t : sheet_tets) {
                Region tet = m_mesh->get<Region>(t);
                std::vector<TCellID> faces = tet.getIDs<Face>();

                for (auto f : faces) {
                    if (linkerG_T.getGeomDim<Face>(f) == 3) {

                        if (linkerG_T.getGeomId<Face>(f) == m_boundarySheets_to_surf[dual_sheets[i].getID()]) {

                            Face face = m_mesh->get<Face>(f);
                            std::vector<TCellID> nodes = face.getIDs<Node>();
                            Node center = m_mesh->newNode(face.center());
                            //std::cout<<"Center node id = "<<center.id()<<std::endl;

                            Region ghost_tet = m_mesh->newTet(nodes[0], nodes[1], nodes[2], center.id());
                            //std::cout<<"ghost tet "<<ghost_tet.id()<<std::endl;

                            m_mesh->mark<Region>(ghost_tet.id(), mark_ghost);

                            Edge e0 = m_mesh->newEdge(nodes[0], center.id());
                            center.add<Edge>(e0.id());
                            ghost_tet.add<Edge>(e0);
                            Edge e1 = m_mesh->newEdge(nodes[1], center.id());
                            center.add<Edge>(e1.id());
                            ghost_tet.add<Edge>(e1);
                            Edge e2 = m_mesh->newEdge(nodes[2], center.id());
                            center.add<Edge>(e2.id());
                            ghost_tet.add<Edge>(e2);

                            //std::cout<<"add edges"<<std::endl;
                            std::vector<Edge> edges = face.get<Edge>();
                            for (auto e : edges) {
                                e.add<Region>(ghost_tet);
                                ghost_tet.add<Edge>(e);
                            }

                            //std::cout<<"add nodes"<<std::endl;
                            std::vector<Node> face_nodes = face.get<Node>();
                            for (auto n : face_nodes) {
                                n.add<Region>(ghost_tet);
                                //ghost_tet.add<Node>(n);
                            }
                            //std::cout<<"add face"<<std::endl;
                            face.add<Region>(ghost_tet);
                            ghost_tet.add<Face>(face);

                            if (((*m_sheet_X)[tet.id()] == -1 || (*m_sheet_X)[tet.id()] == 0 ||
                                 (*m_sheet_X)[tet.id()] == id_sheet) &&
                                ((*m_sheet_Y)[tet.id()] == -1 || (*m_sheet_Y)[tet.id()] == 0 ||
                                 (*m_sheet_Y)[tet.id()] == id_sheet) &&
                                ((*m_sheet_Z)[tet.id()] == -1 || (*m_sheet_Z)[tet.id()] == 0 ||
                                 (*m_sheet_Z)[tet.id()] == id_sheet)) {

                                (*m_block)[ghost_tet.id()] = 0;
                            } else if ((*m_sheet_X)[tet.id()] == id_sheet) {

                                (*m_block)[ghost_tet.id()] = -1;
                                (*m_sheet_X)[ghost_tet.id()] = 0;
                                if ((*m_sheet_Y)[tet.id()] != 0 && (*m_sheet_Z)[tet.id()] != 0) {
                                    (*m_sheet_Y)[ghost_tet.id()] = (*m_sheet_Y)[tet.id()];
                                    (*m_sheet_Z)[ghost_tet.id()] = (*m_sheet_Z)[tet.id()];
                                } else if ((*m_sheet_Y)[tet.id()] != 0) {
                                    (*m_sheet_Y)[ghost_tet.id()] = (*m_sheet_Y)[tet.id()];
                                    (*m_sheet_Z)[ghost_tet.id()] = 0;
                                } else if ((*m_sheet_Z)[tet.id()] != 0) {
                                    (*m_sheet_Y)[ghost_tet.id()] = 0;
                                    (*m_sheet_Z)[ghost_tet.id()] = (*m_sheet_Z)[tet.id()];
                                }

                            } else if ((*m_sheet_Y)[tet.id()] == id_sheet) {

                                (*m_block)[ghost_tet.id()] = -1;
                                (*m_sheet_Y)[ghost_tet.id()] = 0;
                                if ((*m_sheet_X)[tet.id()] != 0 && (*m_sheet_Z)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = (*m_sheet_X)[tet.id()];
                                    (*m_sheet_Z)[ghost_tet.id()] = (*m_sheet_Z)[tet.id()];
                                } else if ((*m_sheet_X)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = (*m_sheet_X)[tet.id()];
                                    (*m_sheet_Z)[ghost_tet.id()] = 0;
                                } else if ((*m_sheet_Z)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = 0;
                                    (*m_sheet_Z)[ghost_tet.id()] = (*m_sheet_Z)[tet.id()];
                                }

                            } else if ((*m_sheet_Z)[tet.id()] == id_sheet) {

                                (*m_block)[ghost_tet.id()] = -1;
                                (*m_sheet_Z)[ghost_tet.id()] = 0;
                                if ((*m_sheet_X)[tet.id()] != 0 && (*m_sheet_Y)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = (*m_sheet_X)[tet.id()];
                                    (*m_sheet_Y)[ghost_tet.id()] = (*m_sheet_Y)[tet.id()];
                                } else if ((*m_sheet_X)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = (*m_sheet_X)[tet.id()];
                                    (*m_sheet_Y)[ghost_tet.id()] = 0;
                                } else if ((*m_sheet_Y)[tet.id()] != 0) {
                                    (*m_sheet_X)[ghost_tet.id()] = 0;
                                    (*m_sheet_Y)[ghost_tet.id()] = (*m_sheet_Y)[tet.id()];
                                }
                            }


                        }
                    }
                }
                //std::cout<<std::endl;
            }
            //}
        }
    }




}
/*----------------------------------------------------------------------------*/
std::vector<int >DualBlockingSession::createCorner(int AZone, std::tuple<int, int, int, int, int, int, int, int> ABlock, std::vector<int> ADual) {

    std::vector<int> block;
    block.push_back(std::get<0>(ABlock));
    block.push_back(std::get<1>(ABlock));
    block.push_back(std::get<2>(ABlock));
    block.push_back(std::get<3>(ABlock));
    block.push_back(std::get<4>(ABlock));
    block.push_back(std::get<5>(ABlock));
    block.push_back(std::get<6>(ABlock));
    block.push_back(std::get<7>(ABlock));

    std::vector<int> corner;
    std::set<int> corner_set;

    //std::cout<<"Corner of zone "<<AZone<<std::endl;
    int X = ADual[0];
    //std::cout<<"X = "<<X<<std::endl;
    int Y = ADual[1];
    //std::cout<<"Y = "<<Y<<std::endl;
    int Z = ADual[2];
    //std::cout<<"Z = "<<Z<<std::endl;

    std::vector<Region> tet_edge;
    std::set<Region> tet_edge_set;

    for (auto r : m_mesh->regions()) {
        if(m_part->value(r) == 1) {
            if ((*m_block)[r] == AZone) {

                Region tet_r = m_mesh->get<Region>(r);
                m_mesh->mark<Region>(r, mark_visited);


                std::vector<Edge> edges = tet_r.get<Edge>();
                for (auto e : edges) {
                    std::vector<Region> temp = e.get<Region>();
                    for (auto t : temp) {
                        tet_edge_set.insert(t);
                    }
                }
                //std::cout<<"nb tet "<<tet_edge_set.size()<<std::endl;
                while (!tet_edge_set.empty()) {

                    tet_edge.insert(tet_edge.end(), tet_edge_set.begin(), tet_edge_set.end());
                    tet_edge_set.clear();

                    for (auto tet_e : tet_edge) {
                        if (m_part->value(tet_e.id()) == 1) {

                            (*m_propagation_round)[tet_e.id()] = 1;

                            //std::cout<<"try for"<<std::endl;


                            if (!m_mesh->isMarked<Region>(tet_e.id(), mark_visited)) {
                                m_mesh->mark<Region>(tet_e.id(), mark_visited);


                                //if more than one sheet
                                if (((*m_sheet_X)[tet_e.id()] != 0 &&
                                     ((*m_sheet_Y)[tet_e.id()] != 0 || (*m_sheet_Z)[tet_e.id()] != 0))
                                    ||
                                    ((*m_sheet_Y)[tet_e.id()] != 0 &&
                                     ((*m_sheet_X)[tet_e.id()] != 0 || (*m_sheet_Z)[tet_e.id()] != 0))
                                    ||
                                    ((*m_sheet_Z)[tet_e.id()] != 0 &&
                                     ((*m_sheet_X)[tet_e.id()] != 0 || (*m_sheet_Y)[tet_e.id()] != 0))) {

                                    //std::cout<<"Tet on multiple sheets"<<std::endl;
                                    //exit(1);
                                    continue;
                                }
                                    //else if another zone within the block
                                else if (std::find(block.begin(), block.end(), (*m_block)[tet_e.id()]) != block.end() &&
                                         (*m_block)[tet_e.id()] != AZone) {
                                    corner_set.insert((*m_block)[tet_e.id()]);
                                    //std::cout<<"Zone "<<(*m_block)[tet_e.id()]<<"found "<<std::endl;
                                }
                                    //else if the tet is in sheet who doesnt belong to the dual of the block
                                else if ((((*m_sheet_X)[tet_e.id()] != X && (*m_sheet_X)[tet_e.id()] != Y &&
                                           (*m_sheet_X)[tet_e.id()] != Z) && (*m_sheet_X)[tet_e.id()] != 0) ||
                                         (((*m_sheet_Y)[tet_e.id()] != X && (*m_sheet_Y)[tet_e.id()] != Y &&
                                           (*m_sheet_Y)[tet_e.id()] != Z) && (*m_sheet_Y)[tet_e.id()] != 0) ||
                                         (((*m_sheet_Z)[tet_e.id()] != X && (*m_sheet_Z)[tet_e.id()] != Y &&
                                           (*m_sheet_Z)[tet_e.id()] != Z) && (*m_sheet_Z)[tet_e.id()] != 0)) {
                                    //exit(2);
                                    continue;
                                }
                                    // if the tet belongs to the same dual zone
                                else if ((*m_block)[tet_e.id()] == AZone) {
                                    //std::cout<<"same zone"<<std::endl;
                                    std::vector<Edge> edge_tet = tet_e.get<Edge>();
                                    for (auto e_t : edge_tet) {
                                        std::vector<Region> temp = e_t.get<Region>();
                                        for (auto t : temp) {
                                            tet_edge_set.insert(t);
                                        }
                                    }
                                } else if ((*m_block)[tet_e.id()] == -1) {
                                    if (m_mesh->isMarked<Region>(tet_e.id(), mark_ghost)) {

                                        bool zone_found = false;
                                        std::vector<Node> nodes_tet = tet_e.get<Node>();
                                        for (auto n : nodes_tet) {
                                            std::vector<Region> tet_nodes = n.get<Region>();
                                            for (auto t_n : tet_nodes) {
                                                if (m_part->value(t_n.id()) == 1) {
                                                    //if(AZone == 2)
                                                    //std::cout<<"test"<<std::endl;
                                                    if (!m_mesh->isMarked<Region>(t_n.id(), mark_visited) &&
                                                        std::find(block.begin(), block.end(), (*m_block)[t_n.id()]) !=
                                                        block.end() &&
                                                        (*m_block)[t_n.id()] != AZone) {

                                                        corner_set.insert((*m_block)[t_n.id()]);
                                                        zone_found = true;
                                                    }
                                                    if (zone_found) {
                                                        break;
                                                    }
                                                }
                                            }
                                            if (zone_found) {
                                                break;
                                            }
                                        }
                                        if (!zone_found) {
                                            std::vector<Edge> edges_tet = tet_e.get<Edge>();
                                            for (auto e_t : edges_tet) {
                                                std::vector<Region> tet_edges = e_t.get<Region>();
                                                for (auto t_e : tet_edges) {
                                                    if (m_part->value(t_e.id()) == 1) {
                                                        if (m_mesh->isMarked<Region>(t_e.id(), mark_ghost)) {
                                                            tet_edge_set.insert(t_e);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        std::vector<Edge> edges_tet = tet_e.get<Edge>();
                                        for (auto e : edges_tet) {
                                            std::vector<Region> tet_edges = e.get<Region>();
                                            for (auto t_e : tet_edges) {
                                                if (m_part->value(t_e.id()) == 1) {
                                                    if (!m_mesh->isMarked<Region>(t_e.id(), mark_visited) &&
                                                        std::find(block.begin(), block.end(), (*m_block)[t_e.id()]) !=
                                                        block.end() &&
                                                        (*m_block)[t_e.id()] != AZone) {
                                                        corner_set.insert((*m_block)[t_e.id()]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                    //else
                                else {
                                    std::cout << "from zone " << AZone << std::endl;
                                    std::cout << "tet error " << tet_e.id() << std::endl;
                                    exit(28);
                                }
                            }
                        }
                    }

                    tet_edge.clear();
                }

                if (corner_set.size() == 3) {
                    for (auto c : corner_set) {
                        corner.push_back(c);
                    }
                    break;
                } else {
                    return corner;
                }
            }
        }
    }

    m_mesh->unmarkAll<Region>(mark_visited);

    return corner;
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::classifyEdges(int AMode) {
    if(AMode == 0) {
        for (auto f : m_mesh->faces()) {
            if (linkerG_T.getGeomDim<Face>(f) == 3) {
                std::vector<Edge> edges = m_mesh->get<Face>(f).get<Edge>();
                for (auto e : edges) {
                    if (linkerG_T.getGeomDim(e.get<Node>()[0]) != 2 && linkerG_T.getGeomDim(e.get<Node>()[1]) != 2) {
                        if (linkerG_T.getGeomId(e.get<Node>()[0]) == linkerG_T.getGeomId(e.get<Node>()[1])) {
                            linkerG_T.linkEdgeToSurface(e.id(), linkerG_T.getGeomId(e.get<Node>()[0]));
                        } else if (linkerG_T.getGeomDim(e.get<Node>()[0]) != 3 &&
                                   linkerG_T.getGeomDim(e.get<Node>()[1]) > 0) {
                            linkerG_T.linkEdgeToSurface(e.id(), linkerG_T.getGeomId(e.get<Node>()[0]));
                        } else {
                            linkerG_T.linkEdgeToSurface(e.id(), linkerG_T.getGeomId(e.get<Node>()[1]));
                        }
                    }
                }
            }
        }
    } else if (AMode > 0){

        for(auto e : m_hmesh->edges()){

            std::vector<int> on_curve_0;
            std::vector<int> on_surf_0;
            std::vector<int> on_volume_0;

            std::vector<int> on_curve_1;
            std::vector<int> on_surf_1;
            std::vector<int> on_volume_1;

            Edge edge = m_hmesh->get<Edge>(e);

            TCellID n0 = edge.getIDs<Node>()[0];
            TCellID n1 = edge.getIDs<Node>()[1];

            int dual_zone_id_0 = 0;
            int dual_zone_id_1 = 0;

            for(auto zone_node : m_dual_zones){
                if(zone_node.second.id() == n0){
                    dual_zone_id_0 = zone_node.first;
                }else if(zone_node.second.id() == n1){
                    dual_zone_id_1 = zone_node.first;
                }
            }

            std::cout<<"Node id "<<n0<<" dual zone "<<dual_zone_id_0<<std::endl;
            std::cout<<"Node id "<<n1<<" dual zone "<<dual_zone_id_1<<std::endl;

            if(dual_zone_id_0 == 0 || dual_zone_id_1 == 0) throw GMDSException("ERROR: In edge block classification : "
                                                                               "no dual zone corresponding to block node");

            for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_0]){
                if(geom_infos.first == 2){
                    on_curve_0.push_back(geom_infos.second);
                }else if(geom_infos.first == 3){
                    on_surf_0.push_back(geom_infos.second);
                }else{
                    on_volume_0.push_back(geom_infos.second);
                }
            }
            for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_1]){
                if(geom_infos.first == 2){
                    on_curve_1.push_back(geom_infos.second);
                }else if(geom_infos.first == 3){
                    on_surf_1.push_back(geom_infos.second);
                }else{
                    on_volume_1.push_back(geom_infos.second);
                }
            }

            int geom_dim = 0;
            int geom_id = 0;
            bool common_curve_found = false;

            if(!on_curve_0.empty() && !on_curve_1.empty()){
                for(auto c0_id : on_curve_0){
                    for(auto c1_id : on_curve_1){
                        if(c0_id == c1_id){
                            geom_dim = 2;
                            geom_id = c0_id;
                            common_curve_found = true;
                        }
                    }
                }
            }
            int cpt_surfs = 0;
            if(!common_curve_found && (!on_surf_0.empty() && !on_surf_1.empty())){
                for(auto s0_id : on_surf_0){
                    for(auto s1_id : on_surf_1){
                        if(s0_id == s1_id){
                            geom_dim = 3;
                            geom_id = s0_id;
                            cpt_surfs++;
                        }
                    }
                }
            }
            if(cpt_surfs == 2 && !common_curve_found){
                geom_dim = 0;
                geom_id = 1;
                std::cout<<"Dual zone "<<dual_zone_id_0<<std::endl;
                for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_0]){
                    std::cout<<"Dim = "<<geom_infos.first<<", Id = "<<geom_infos.second<<std::endl;
                }
                std::cout<<"Dual zone "<<dual_zone_id_1<<std::endl;
                for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_1]){
                    std::cout<<"Dim = "<<geom_infos.first<<", Id = "<<geom_infos.second<<std::endl;
                }

                /*throw GMDSException("ERROR: In edge block classification : "
                                    "two surfaces in common but no curve found in "+std::to_string(dual_zone_id_0)+
                                    " and "+std::to_string(dual_zone_id_1));*/
            }

            if(geom_dim == 2){
                linkerH_G.linkToCurve(edge, geom_id);
            }else if(geom_dim == 3){
                linkerH_G.linkToSurface(edge, geom_id);
            }else{
                //link to volume but for now it equal to no link
            }

        }



/*
        for(auto e : m_hmesh->edges()){
            Edge edge = m_hmesh->get<Edge>(e);
            int geom_ID_n0 = linkerH_G.getGeomId<Node>(edge.getIDs<Node>()[0]);
            int geom_Dim_n0 = linkerH_G.getGeomDim<Node>(edge.getIDs<Node>()[0]);

            int geom_ID_n1 = linkerH_G.getGeomId<Node>(edge.getIDs<Node>()[1]);
            int geom_Dim_n1 = linkerH_G.getGeomDim<Node>(edge.getIDs<Node>()[1]);

            if(geom_Dim_n0 == cad::GeomMeshLinker::LINK_POINT){
                cad::GeomPoint* p0 = manager.getPoint(geom_ID_n0);
                if(geom_Dim_n1 == cad::GeomMeshLinker::LINK_POINT){
                    cad::GeomPoint* p1 = manager.getPoint(geom_ID_n1);
                    //NN
                    int common_edge = manager.getCommonCurve(p0,p1);
                    int common_face = manager.getCommonSurface(p0,p1);
                    int cpt_common_faces = 0;
                    for(auto s0 : p0->surfaces()){
                        for(auto s1 : p1->surfaces()) {
                            if (s0->id() == s1->id()) {
                                common_face = s0->id();
                                cpt_common_faces++;
                            }
                        }
                    }
                    if(common_edge != -1){
                        //NN -> E
                        //On a trouvé la courbe commune entre les deux points donc l'arete de bloc est sur une courbe
                        linkerH_G.linkToCurve(edge,common_edge);
                    } else if(common_face != -1 && cpt_common_faces == 1){
                        //NN -> F
                        //On a trouvé la surface commune entre les deux points donc l'arete de bloc est sur une courbe
                        linkerH_G.linkToSurface(edge,common_face);
                    }


                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_CURVE){
                    //NE
                    bool edge_found = false;
                    //On cherche si N0 est une extremité de C1
                    for(auto c : p0->curves()){
                        if(c->id() == geom_ID_n1){
                            //NE -> E
                            linkerH_G.linkToCurve(edge,geom_ID_n1);
                            edge_found = true;
                        }
                    }
                    if(!edge_found){
                        //N0 n'est pas une extremité de C1 on cherche alors si les deux on une surface commune
                        cad::GeomCurve* c1 = manager.getCurve(geom_ID_n1);
                        for(auto s0 : p0->surfaces()){
                            for(auto s1 : c1->surfaces()) {
                                if (s0->id() == s1->id()) {
                                    //NE -> F
                                    linkerH_G.linkToSurface(edge, s1->id());
                                }
                            }
                        }
                    }

                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_SURFACE){
                    //NF
                    for(auto s : p0->surfaces()){
                        if(s->id() == geom_ID_n1){
                            //NF -> F
                            linkerH_G.linkToSurface(edge,geom_ID_n1);
                        }
                    }
                }
            }else if(geom_Dim_n0 == cad::GeomMeshLinker::LINK_CURVE){
                cad::GeomCurve* c0 = manager.getCurve(geom_ID_n0);
                if(geom_Dim_n1 == cad::GeomMeshLinker::LINK_POINT){
                    cad::GeomPoint* p1 = manager.getPoint(geom_ID_n1);
                    //EN
                    bool edge_found = false;
                    //On cherche si N1 est une extremité de C0
                    for(auto c : p1->curves()){
                        if(c->id() == geom_ID_n0){
                            //EN -> E
                            linkerH_G.linkToCurve(edge,geom_ID_n0);
                            edge_found = true;
                        }
                    }
                    if(!edge_found){
                        //N1 n'est pas une extremité de C0 on cherche alors si les deux on une surface commune
                        for(auto s1 : p1->surfaces()){
                            for(auto s0 : c0->surfaces()) {
                                if (s1->id() == s0->id()) {
                                    //EN -> F
                                    linkerH_G.linkToSurface(edge, s0->id());
                                }
                            }
                        }
                    }


                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_CURVE){
                    //EE
                    cad::GeomCurve* c1 = manager.getCurve(geom_ID_n1);
                    if(geom_ID_n0 == geom_ID_n1){
                        linkerH_G.linkToCurve(edge,geom_ID_n0);
                    } else{
                        int common_face = manager.getCommonSurface(c0,c1);
                        int common_face_test = 0;
                        for(auto s0 : c0->surfaces()){
                            for(auto s1 : c1->surfaces()){
                                if(s0->id() == s1->id()){
                                    common_face_test = s0->id();
                                }
                            }
                        }
                        //if(common_face != common_face_test) throw GMDSException("ERROR SUR LA SURFACE COMMUNE ENTRE DEUX COURBES");
                        if(common_face_test != -1){
                            linkerH_G.linkToSurface(edge, common_face_test);
                        }
                    }

                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_SURFACE){
                    //EF
                    for(auto s : c0->surfaces()){
                        if(s->id() == geom_ID_n1){
                            //EF -> F
                            linkerH_G.linkToSurface(edge,geom_ID_n1);
                        }
                    }
                }
            }else if(geom_Dim_n0 == cad::GeomMeshLinker::LINK_SURFACE){
                cad::GeomSurface* s0 = manager.getSurface(geom_ID_n1);
                if(geom_Dim_n1 == cad::GeomMeshLinker::LINK_POINT){
                    //FN
                    cad::GeomPoint* p1 = manager.getPoint(geom_ID_n1);

                    for(auto s : p1->surfaces()){
                        if(s->id() == geom_ID_n0){
                            //FN -> F
                            linkerH_G.linkToSurface(edge,geom_ID_n0);
                        }
                    }
                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_CURVE){
                    //FE
                    cad::GeomCurve* c1 = manager.getCurve(geom_ID_n1);
                    for(auto s : c1->surfaces()){
                        if(s->id() == geom_ID_n0){
                            //FE -> F
                            linkerH_G.linkToSurface(edge,geom_ID_n0);
                        }
                    }

                } else if (geom_Dim_n1 == cad::GeomMeshLinker::LINK_SURFACE){
                    //FF
                    if(geom_ID_n0 == geom_ID_n1){
                        linkerH_G.linkToSurface(edge,geom_ID_n0);
                    }
                }
            }
        }*/
    }
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::classifyBlockFaces(){

    for(auto f : m_hmesh->faces()){

        std::vector<int> on_curve_0;
        std::vector<int> on_surf_0;
        std::vector<int> on_volume_0;

        std::vector<int> on_curve_1;
        std::vector<int> on_surf_1;
        std::vector<int> on_volume_1;

        std::vector<int> on_curve_2;
        std::vector<int> on_surf_2;
        std::vector<int> on_volume_2;

        std::vector<int> on_curve_3;
        std::vector<int> on_surf_3;
        std::vector<int> on_volume_3;

        Face face = m_hmesh->get<Face>(f);

        TCellID n0 = face.getIDs<Node>()[0];
        TCellID n1 = face.getIDs<Node>()[1];
        TCellID n2 = face.getIDs<Node>()[2];
        TCellID n3 = face.getIDs<Node>()[3];

        int dual_zone_id_0 = 0;
        int dual_zone_id_1 = 0;
        int dual_zone_id_2 = 0;
        int dual_zone_id_3 = 0;

        for(auto zone_node : m_dual_zones){
            if(zone_node.second.id() == n0){
                dual_zone_id_0 = zone_node.first;
            }else if(zone_node.second.id() == n1){
                dual_zone_id_1 = zone_node.first;
            }else if(zone_node.second.id() == n2){
                dual_zone_id_2 = zone_node.first;
            }else if(zone_node.second.id() == n3){
                dual_zone_id_3 = zone_node.first;
            }
        }

        if(dual_zone_id_0 == 0 || dual_zone_id_1 == 0 || dual_zone_id_2 == 0 || dual_zone_id_3 == 0) throw GMDSException("ERROR: In edge block classification : "
                                                                                                                         "no dual zone corresponding to block node");

        for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_0]){
            if(geom_infos.first == 2){
                on_curve_0.push_back(geom_infos.second);
            }else if(geom_infos.first == 3){
                on_surf_0.push_back(geom_infos.second);
            }else{
                on_volume_0.push_back(geom_infos.second);
            }
        }
        for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_1]){
            if(geom_infos.first == 2){
                on_curve_1.push_back(geom_infos.second);
            }else if(geom_infos.first == 3){
                on_surf_1.push_back(geom_infos.second);
            }else{
                on_volume_1.push_back(geom_infos.second);
            }
        }
        for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_2]){
            if(geom_infos.first == 2){
                on_curve_2.push_back(geom_infos.second);
            }else if(geom_infos.first == 3){
                on_surf_2.push_back(geom_infos.second);
            }else{
                on_volume_2.push_back(geom_infos.second);
            }
        }
        for(auto geom_infos : m_dual_zones_classifications[dual_zone_id_3]){
            if(geom_infos.first == 2){
                on_curve_3.push_back(geom_infos.second);
            }else if(geom_infos.first == 3){
                on_surf_3.push_back(geom_infos.second);
            }else{
                on_volume_3.push_back(geom_infos.second);
            }
        }

        int geom_dim = 0;
        int geom_id = 0;



        if(!on_surf_0.empty() && !on_surf_1.empty() && !on_surf_2.empty() && !on_surf_3.empty()){
            for(auto s0_id : on_surf_0){
                for(auto s1_id : on_surf_1){
                    for(auto s2_id : on_surf_2){
                        for(auto s3_id : on_surf_3){
                            if(s0_id == s1_id && s1_id == s2_id && s2_id == s3_id && s3_id == s0_id){
                                geom_dim = 3;
                                geom_id = s0_id;
                            }
                        }
                    }
                }
            }
        }

        std::vector<int> e_on_curve;
        bool is_a_curve_conv = false;
        if(geom_dim == 3){
            // meme si les 4 sommets on une surface en commun, on peut quand meme etre une face volumique
            for(auto e : face.getIDs<Edge>()){
                if(linkerH_G.getGeomDim<Edge>(e) == 2){
                    e_on_curve.push_back(linkerH_G.getGeomId<Edge>(e));
                }
            }
            if(e_on_curve.size() >= 3){
                for(auto c_id : e_on_curve){
                    cad::GeomCurve* c = manager.getCurve(c_id);
                    if(c->getCurvatureInfo() == cad::GeomCurve::CONVEX){
                        is_a_curve_conv = true;
                        break;
                    }
                }
                //Si une des courbes incidente est convexe alors on ne peut pas être volumique
                if(!is_a_curve_conv){
                    geom_dim = 0;
                    geom_id = 0;
                }
            }
        }

        std::cout<<"Face "<<f<<std::endl;
        std::cout<<"Dim = "<<geom_dim<<std::endl;
        std::cout<<"Id = "<<geom_id<<std::endl;
        std::cout<<"Surfs from dual region "+std::to_string(dual_zone_id_0)<<std::endl;
        for(auto s : on_surf_0){
            std::cout<<s<<",";
        }
        std::cout<<std::endl;
        std::cout<<"Surfs from dual region "+std::to_string(dual_zone_id_1)<<std::endl;
        for(auto s : on_surf_1){
            std::cout<<s<<",";
        }
        std::cout<<std::endl;
        std::cout<<"Surfs from dual region "+std::to_string(dual_zone_id_2)<<std::endl;
        for(auto s : on_surf_2){
            std::cout<<s<<",";
        }
        std::cout<<std::endl;
        std::cout<<"Surfs from dual region "+std::to_string(dual_zone_id_3)<<std::endl;
        for(auto s : on_surf_3){
            std::cout<<s<<",";
        }
        std::cout<<std::endl;
        std::cout<<"Nb edges of face on curve "<<e_on_curve.size()<<std::endl;
        if(!e_on_curve.empty()){
            for(auto e : e_on_curve){
                std::cout<<e<<",";
            }
            std::cout<<std::endl;
            if(is_a_curve_conv)
                std::cout<<"a curve is convex"<<std::endl;
        }


        if(geom_dim == 3){
            linkerH_G.linkToSurface(face, geom_id);
        }

    }

    /*
    for(auto f : m_hmesh->faces()){
        Face face = m_hmesh->get<Face>(f);
        std::vector<TCellID> edge_f = face.getIDs<Edge>();
        std::vector<int> e_on_curve;
        std::vector<int> e_on_surface;
        std::vector<int> e_on_volume;
        for(auto e : edge_f){
            auto info = linkerH_G.getGeomInfo<Edge>(e);
            if(info.first == 2){
                std::cout<<"On curve "<<info.second<<std::endl;
                e_on_curve.push_back(info.second);
            } else if(info.first == 3) {
                std::cout<<"On surf "<<info.second<<std::endl;
                e_on_surface.push_back(info.second);
            } else{
                std::cout<<"In volume"<<std::endl;
                e_on_volume.push_back(0);
            }
        }
        if(e_on_curve.size()+ e_on_surface.size()+e_on_volume.size() != 4) throw GMDSException("Wrong number of edge classifications");
        if(e_on_volume.empty()){
           if(e_on_surface.size() == 4){
               //FFFF
               bool s_0_1 = e_on_surface[0] == e_on_surface[1];
               bool s_1_2 = e_on_surface[1] == e_on_surface[2];
               bool s_2_3 = e_on_surface[2] == e_on_surface[3];
               if(s_0_1 && s_1_2 && s_2_3){
                   linkerH_G.linkToSurface(face,e_on_surface[0]);
               }
           } else if (e_on_surface.size() == 3){
               //FFFE
               bool s_0_1 = e_on_surface[0] == e_on_surface[1];
               bool s_1_2 = e_on_surface[1] == e_on_surface[2];
               if(s_0_1 && s_1_2){
                   linkerH_G.linkToSurface(face,e_on_surface[0]);
               }
           } else if (e_on_surface.size() == 2){
               //FFEE
               if(e_on_surface[0] == e_on_surface[1]){
                   linkerH_G.linkToSurface(face,e_on_surface[0]);
               }
           } else if(e_on_surface.size() == 1){
               //FEEE
               cad::GeomCurve* c0 = manager.getCurve(e_on_curve[0]);
               cad::GeomCurve* c1 = manager.getCurve(e_on_curve[1]);
               cad::GeomCurve* c2 = manager.getCurve(e_on_curve[2]);

               std::vector<cad::GeomCurve*> curves;
               curves.push_back(c0);
               curves.push_back(c1);
               curves.push_back(c2);

               cad::GeomSurface* s = manager.getSurface(e_on_surface[0]);

               int cpt_adjacent = 0;
               for(auto c : curves){
                   if(std::find(c->surfaces().begin(), c->surfaces().end(), s) != c->surfaces().end()){
                       cpt_adjacent++;
                   }
               }
               if(cpt_adjacent == 3){
                   linkerH_G.linkToSurface(face,e_on_surface[0]);
               }

           } else if(e_on_curve.size() == 4){
               //EEEE
               cad::GeomCurve* c0 = manager.getCurve(e_on_curve[0]);
               cad::GeomCurve* c1 = manager.getCurve(e_on_curve[1]);
               cad::GeomCurve* c2 = manager.getCurve(e_on_curve[2]);
               cad::GeomCurve* c3 = manager.getCurve(e_on_curve[3]);

               int common_surf_01 = manager.getCommonSurface(c0,c1);
               int common_surf_02 = manager.getCommonSurface(c0,c2);
               int common_surf_03 = manager.getCommonSurface(c0,c3);

               bool common_012 = common_surf_01 == common_surf_02;
               bool common_013 = common_surf_01 == common_surf_03;
               if(common_012 && common_013){
                   linkerH_G.linkToSurface(face,common_surf_01);
               }
           } else{
               throw GMDSException("ERROR Classification face");
           }
           std::cout<<"Face "<<f<<" on surface "<<linkerH_G.getGeomId<Face>(f)<<std::endl;
        }else{
            std::cout<<"Face "<<f<<" in volume"<<std::endl;
        }




        //======================================================
        std::set<int> surface_ids;
        for(auto n : face.getIDs<Node>()){
            if(linkerH_G.getGeomDim<Node>(n) == 3){
                surface_ids.insert(linkerH_G.getGeomId<Node>(n));
            }
        }
        if(surface_ids.size() == 1){
            int id = *surface_ids.begin();
            linkerH_G.linkToSurface(face,id);
        }
    }*/
}
/*----------------------------------------------------------------------------*/
bool DualBlockingSession::deleteSheet(int AID){
    int i =0;
    if(dual_sheets.size()==1){
        for (auto t:dual_sheets[0].getSurface()) {
            uncolorTet(t, AID);
            if (!isColored(t)) { (*m_block)[t] = 0; }
        }
        dual_sheets.clear();
    }else {
        for (auto s:dual_sheets) {
            if (s.getID() == AID) {
                for (auto t:s.getSurface()) {
                    uncolorTet(t, AID);
                    if (!isColored(t)) { (*m_block)[t] = 0; }
                }
                dual_sheets.erase(dual_sheets.begin() + i  );
            }
            i++;
        }
    }
}
/*----------------------------------------------------------------------------*/
int DualBlockingSession::getSheetOfTet(TCellID AID){
    if((*m_sheet_X)[AID] > 0) {return (*m_sheet_X)[AID];}
    else if((*m_sheet_Y)[AID] > 0) {return (*m_sheet_Y)[AID];}
    else if((*m_sheet_Z)[AID] > 0) {return (*m_sheet_Z)[AID];}
    else{return 0;}
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::uncolorTet(gmds::TCellID ATetID, int ASheetID){

    if((*m_sheet_X)[ATetID] == ASheetID) {(*m_sheet_X)[ATetID] = 0;}
    if((*m_sheet_Y)[ATetID] == ASheetID) {(*m_sheet_Y)[ATetID] = 0;}
    if((*m_sheet_Z)[ATetID] == ASheetID) {(*m_sheet_Z)[ATetID] = 0;}
}
/*----------------------------------------------------------------------------*/
bool DualBlockingSession::getFrameAxis(double result[][3],TCellID AID){

    //We suppose for now that we dont select a cell inside the volume

    math::Point center;
    Region r =  m_mesh->get<Region>(AID);
    math::Vector3d norm;
    bool found = false;
    for(auto f : r.get<Face>()){
        std::vector<TCellID> nodes = f.getIDs<Node>();
        if(linkerG_T.getGeomDim<Node>(nodes[0]) == 3 && linkerG_T.getGeomDim<Node>(nodes[1]) == 3 && linkerG_T.getGeomDim<Node>(nodes[2]) == 3){
            std::cout<<"Face au bord"<<std::endl;
            center = f.center();
            norm = f.normal();
            found = true;
            break;
        }
    }
    if(!found) return false;

    double max_length = 0;
    for(auto e : m_mesh->edges()){
        Edge edge = m_mesh->get<Edge>(e);
        if(max_length < edge.length()){
            max_length = edge.length();
        }
    }
    int x = 0;
    int y = 0;
    int i;
    for(i = 0; i<3; i++){
        std::cout<<"Dot = "<<fabs((*m_tetra_chart)[AID][i].dot(norm))<<std::endl;
        if(fabs((*m_tetra_chart)[AID][i].dot(norm)) > 0.9){
            break;
        }
    }
    std::cout<<"i = "<<i<<std::endl;
    if(i == 0){
        x = 1; y = 2;
    }else if(i == 1){
        x = 0; y = 2;
    }else{
        x = 0; y = 1;
    }
    //We store the indices of x and y to use it later for the surface creation normal
    m_axis[0] = x;
    m_axis[1] = y;

    result[0][0] = center.X();
    result[0][1] = center.Y();
    result[0][2] = center.Z();

    result[1][0] = center.X()+((*m_tetra_chart)[AID][x].X()*(max_length));
    result[1][1] = center.Y()+((*m_tetra_chart)[AID][x].Y()*(max_length));
    result[1][2] = center.Z()+((*m_tetra_chart)[AID][x].Z()*(max_length));

    result[2][0] = center.X()+((*m_tetra_chart)[AID][y].X()*(max_length));
    result[2][1] = center.Y()+((*m_tetra_chart)[AID][y].Y()*(max_length));
    result[2][2] = center.Z()+((*m_tetra_chart)[AID][y].Z()*(max_length));

    return true;

}
/*----------------------------------------------------------------------------*/
Mesh* DualBlockingSession::getBlocks() {
    return m_hmesh;
}
/*----------------------------------------------------------------------------*/
gmds::Mesh* DualBlockingSession::getSurfaceMesh(int ASheetID) {
    for(auto ds : dual_sheets){
        if(ds.getID() == ASheetID){
            return ds.getSurfaceMesh();
        }
    }
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::resetDual() {
    for(auto r : m_mesh->regions()){
        if(m_block->value(r) != -1){
            m_block->set(r,0);
        }
    }
    m_zone_border.clear();
    m_dual_zones.clear();
    m_hmesh->clear();
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::hardResetDual() {
    for(auto r : m_mesh->regions()){
        m_sheet_X->set(r,0);
        m_sheet_Y->set(r,0);
        m_sheet_Z->set(r,0);

        m_block->set(r,0);
    }

    m_zone_border.clear();
    m_dual_zones.clear();
    m_hmesh->clear();
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::resetBlocks() {

    for(auto r : m_hmesh->regions()){
        m_hmesh->deleteRegion(r);
    }
    for(auto f : m_hmesh->faces()){
        m_hmesh->deleteFace(f);
    }
    for(auto e : m_hmesh->edges()){
        m_hmesh->deleteEdge(e);
    }
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> DualBlockingSession::getSingularGraph(){

    std::vector<TCellID> sGraph;
    for(auto r : m_mesh->regions()){
        if((*m_sing)[r] > 0){
            sGraph.push_back(r);
        }
    }
    return sGraph;
}
/*----------------------------------------------------------------------------*/
bool DualBlockingSession::checkDualRegionsValidity(){

    std::map<int,TCellID> bad_dual_regions;

    //Version bourrin nulle mais qui marche

    std::map<int, std::vector<TCellID>> dual_to_tets;
    for(int i = 1; i< m_dual_zones.size()+1; i++){
        std::vector<TCellID> tets;
        dual_to_tets[i] = tets;
    }

    for(auto r : m_mesh->regions()){
        if(m_block->value(r) > 0){
            dual_to_tets[m_block->value(r)].push_back(r);
        }
    }
    for(auto const &region : dual_to_tets){
        bool geom_node_found = false;
        std::set<int> geom_node_ids;
        for(auto r : region.second){
            Region tet = m_mesh->get<Region>(r);
            //Ici on teste si le tétra est une singularité
            if(m_sing->value(r) == 3){
                geom_node_found = true;
            }
            // Ici on teste si on a un tetra au bord et si oui, s'il est singulier
            for(auto const &f : tet.get<Face>()){
                if(f.get<Region>().size() == 1 && m_sing->value(r) == 1){
                    geom_node_found = true;
                }
            }
            // Ici on teste si un des sommets du tetra est un sommet geom
            for(auto const &n : tet.get<Node>()){
                if(linkerG_T.getGeomDim(n) == 1){
                    geom_node_ids.insert(linkerG_T.getGeomId(n));
                    geom_node_found = true;
                }
            }
            //if(geom_node_found) break;
        }

        //La zone n'est pas valide si on n'a aucun sommet géom ou singulier, ou bien si on a plusieurs sommets géom
        if(!geom_node_found || geom_node_ids.size() > 1)
            bad_dual_regions.emplace(region.first,region.second[0]);
    }

    m_wrong_dual_regions = bad_dual_regions;
    colorWrongDualRegions();
    return bad_dual_regions.empty();
}
/*----------------------------------------------------------------------------*/
int DualBlockingSession::removeWrongDualRegions(){
    for(auto r : m_mesh->regions()){
        if(m_block->value(r) == -2){
            m_block->set(r, -1);
        }
    }
    for(auto w : m_wrong_dual_regions){
        m_hmesh->deleteNode(m_dual_zones[w.first]);
        m_dual_zones.erase(w.first);
    }

    return m_dual_zones.size();
}
/*----------------------------------------------------------------------------*/
int DualBlockingSession::recolorDual(){
    for(auto w : m_wrong_dual_regions){
        int dual_region_id = w.first;
        Region origin = m_mesh->get<Region>(w.second);
        m_block->set(w.second, dual_region_id);

        std::vector<Region> regions;

        std::vector<Edge> edges = origin.get<Edge>();
        for(auto e : edges){
            std::vector<Region> e_regions = e.get<Region>();
            for(auto e_r : e_regions){
                if(m_block->value(e_r.id()) == -2){
                    m_block->set(e_r.id(), dual_region_id);
                    regions.push_back(e_r);
                }
            }
        }

        while(!regions.empty()){
            std::vector<Region> tmp;
            for(auto r : regions) {
                std::vector<Edge> edges_r = r.get<Edge>();
                for (auto e : edges_r) {
                    std::vector<Region> e_regions = e.get<Region>();
                    for (auto e_r : e_regions) {
                        if (m_block->value(e_r.id()) == -2) {
                            m_block->set(e_r.id(), dual_region_id);
                            tmp.push_back(e_r);
                        }
                    }
                }
            }
            regions = tmp;
        }
    }

    IGMeshIOService ioService_m(m_mesh);
    VTKWriter vtkWriter_m(&ioService_m);
    vtkWriter_m.setCellOptions(gmds::N|gmds::R);
    vtkWriter_m.setDataOptions(gmds::N|gmds::R);
    vtkWriter_m.write("/home/simon/Data/Results_debug/blocks.vtk");

    return m_dual_zones.size();
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::colorWrongDualRegions(){
    std::vector<int> ids;
    for(auto w : m_wrong_dual_regions){
        ids.push_back(w.first);
    }

    for(auto r : m_mesh->regions()){
        if(std::find(ids.begin(), ids.end(), m_block->value(r)) != ids.end()){
            m_block->set(r, -2);
        }
    }
}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::correctionBlockEdgesClassification(){

    // On récupère les surfaces geom qui ont eu une correction de frame
    Variable<int>* free_bound = m_mesh->getVariable<int,GMDS_NODE>("free_bnd_nodes");
    std::set<int> frame_corrected_surfaces;
    for(auto n : m_mesh->nodes()){
        if(free_bound->value(n) == 1){
            frame_corrected_surfaces.insert(linkerG_T.getGeomId<Node>(n));
        }
    }

    for(auto e : m_hmesh->edges()){
        // On récupère les arêtes de blocs qui sont classifiées sur ces surfaces
        if(linkerH_G.getGeomDim<Edge>(e) == 3 &&
                frame_corrected_surfaces.find(linkerH_G.getGeomId<Edge>(e))
                != frame_corrected_surfaces.end()){

            Edge h_edge = m_hmesh->get<Edge>(e);
            bool f1_found = false;
            Face f1;
            Face f2;

            // On récupère les deux faces au bord adjacentes
            for(auto f : h_edge.get<Face>()){
                if(f.getIDs<Region>().size() == 1){
                    if(!f1_found){
                        f1 = f;
                        f1_found = true;
                    }
                    else{
                        f2 = f;
                    }
                }
            }
            std::cout<<"F1 "<<f1.id()<<std::endl;
            std::cout<<"F2 "<<f2.id()<<std::endl;

            std::vector<Node> ns1 = f1.get<Node>();
            math::Vector3d v1 = f1.normal();
            math::Vector3d v2 = f2.normal();
            math::Point c2  = f2.center();
            bool convex=true;
            if(math::Orientation::orient3d(c2,ns1[0].getPoint(),
                                           ns1[1].getPoint(),
                                           ns1[2].getPoint())==math::Orientation::NEGATIVE){
                //its non-convex
                convex=false;
            }
            double a = math::Constants::INVPIDIV180*v1.angle(v2);
            if(!convex)
                a+=180;
            if(a<180)
                std::cout<<"convex"<<std::endl;
            else if( a> 180)
                std::cout<<"concave"<<std::endl;
        }
    }

}
/*----------------------------------------------------------------------------*/
void DualBlockingSession::testLissageSheet(){
    /*int mark_node_done = m_mesh->newMark<Node>();
    for(auto d : dual_sheets) {
        if (!d.isBoundary()) {
            int sheet_id = d.getID();
            std::vector<TCellID> tets = d.getSurface();
            for (int i = 0; i < 1; i++) {
                for (auto t : tets) {
                    Region r = m_mesh->get<Region>(t);
                    for (auto n : r.get<Node>()) {
                        if(linkerG_T.getGeomDim(n) != 2){

                            if (!m_mesh->isMarked<Node>(n.id(), mark_node_done)) {
                                std::set<Node> neighbors;
                                for (auto f : n.get<Face>()) {
                                    if (f.nbRegions() == 2) {
                                        bool r_f0 = m_sheet_X->value(f.getIDs<Region>()[0]) == sheet_id ||
                                                    m_sheet_Y->value(f.getIDs<Region>()[0]) == sheet_id ||
                                                    m_sheet_Z->value(f.getIDs<Region>()[0]) == sheet_id;
                                        bool r_f1 = m_sheet_X->value(f.getIDs<Region>()[1]) == sheet_id ||
                                                    m_sheet_Y->value(f.getIDs<Region>()[1]) == sheet_id ||
                                                    m_sheet_Z->value(f.getIDs<Region>()[1]) == sheet_id;

                                        if (r_f0 ^ r_f1) {
                                            std::cout<<r_f0<<","<<r_f1<<std::endl;
                                            for (auto n_f : f.get<Node>()) {
                                                if (n_f != n) {
                                                    if(linkerG_T.getGeomDim(n) == 3){
                                                        if(linkerG_T.getGeomDim(n_f) == 3 && linkerG_T.getGeomId(n_f) == linkerG_T.getGeomId(n)){
                                                            neighbors.insert(n_f);
                                                        }
                                                    }else{
                                                        neighbors.insert(n_f);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                std::vector<math::Point> points;
                                points.reserve(neighbors.size());
                                for (auto node : neighbors) {
                                    points.push_back(node.getPoint());
                                }
                                math::Point center = math::Point::massCenter(points);

                                math::Vector3d norm;
                                if (m_sheet_X->value(t) == sheet_id) {
                                    norm = m_tetra_chart->value(t)[0];
                                } else if (m_sheet_Y->value(t) == sheet_id) {
                                    norm = m_tetra_chart->value(t)[1];
                                } else if (m_sheet_Z->value(t) == sheet_id) {
                                    norm = m_tetra_chart->value(t)[2];
                                }

                                math::Plane plane(n.getPoint(), norm);
                                math::Point project = plane.project(center);
                                math::Vector3d vec(project, center);

                                math::Point new_point(n.getPoint().X() + vec.X(), n.getPoint().Y() + vec.Y(),
                                                      n.getPoint().Z() + vec.Z());

                                m_mesh->mark<Node>(n.id(),mark_node_done);
                                if(linkerG_T.getGeomDim(n) == 2){
                                    manager.getCurve(linkerG_T.getGeomId(n))->project(center);
                                } else if (linkerG_T.getGeomDim(n) == 3){
                                    std::cout<<"tet "<<t<<std::endl;
                                    std::cout<<"node "<<n.id()<<std::endl;
                                    for(auto node : neighbors){
                                        std::cout<<node.id()<<", ";
                                    }
                                    std::cout<<std::endl;
                                    manager.getSurface(linkerG_T.getGeomId(n))->project(center);
                                }
                                //n.setPoint(center);
                            }
                        }
                    }
                }
            }
            m_mesh->unmarkAll<Node>(mark_node_done);
        }
    }
    m_mesh->unmarkAll<Node>(mark_node_done);
    m_mesh->freeMark<Node>(mark_node_done);*/
}