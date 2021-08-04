/*----------------------------------------------------------------------------*/
#include <medusa/view/GraphicView.h>
/*----------------------------------------------------------------------------*/
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkCellData.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <gmds/utils/Exception.h>
#include <medusa/model/MedusaBackEnd.h>
#include <medusa/control/MouseInteractorStyle.h>
#include <vtkThreshold.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkLine.h>
#include <vtkSphereSource.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
GraphicView::GraphicView(const std::string& AName, GraphicView::ViewType AType)
        : View(), m_name(AName), m_type(AType)
{
    m_visible=false;

    //We initialize a picking listener ont this view
    m_vtk_renderer = vtkRenderer::New();
    m_vtk_renderer->SetBackground( 0.92, 0.92, 0.92 );

    m_vtk_window = vtkRenderWindow::New();
    m_vtk_window->SetSize( 1000, 800 );
    m_vtk_window->AddRenderer( m_vtk_renderer );
    m_vtk_window->SetWindowName(m_name.c_str());
    m_vtk_interactor= vtkRenderWindowInteractor::New();
    m_vtk_interactor->SetRenderWindow(m_vtk_window);
    //Initialization is mandatory to be displayed
    m_vtk_interactor->Initialize();

    // Set the custom stype to use for interaction.
    m_mouse=
            vtkSmartPointer<MouseInteractorStyle>::New();
    m_mouse->SetDefaultRenderer(m_vtk_renderer);
    m_vtk_interactor->SetInteractorStyle(m_mouse);

    selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    selectedActor  = vtkSmartPointer<vtkActor>::New();

    m_vtk_cut_actor = vtkSmartPointer<vtkActor>::New();

    m_vtk_tets_set = vtkSmartPointer<vtkActor>::New();
    m_vtk_renderer->AddActor(m_vtk_tets_set);

    m_mapper_dual = vtkDataSetMapper::New();
    m_mapper_dual->CreateDefaultLookupTable();

    m_mapper_dual->GetLookupTable()->Build();

    m_textActor = vtkSmartPointer<vtkTextActor>::New();
    m_textActor->SetPosition2 ( 10, 40 );
    m_textActor->GetTextProperty()->SetFontSize ( 24 );
    m_textActor->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );
    m_vtk_renderer->AddActor2D ( m_textActor );

}
/*----------------------------------------------------------------------------*/
GraphicView::~GraphicView() {
    m_vtk_renderer->Delete();
    m_vtk_window->Delete();
    m_vtk_interactor->Delete();
   for(auto a: m_vtk_actors)
            a->Delete();

}
/*----------------------------------------------------------------------------*/
void GraphicView::setVisibleOFF() {m_visible=false;}
/*----------------------------------------------------------------------------*/
void GraphicView::setVisibleON() {m_visible=true;}
/*----------------------------------------------------------------------------*/
void GraphicView::refresh() {

    setVisibleON();
    //Render the scene .
    m_vtk_window->Render();

}
/*----------------------------------------------------------------------------*/
void GraphicView::selectVTKCell(int ADIm, vtkIdType ACellID) {
    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();


    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    ids->InsertNextValue(ACellID);

    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);

    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());


    selectedMapper->SetInputData(selected);
    selectedActor->SetMapper(selectedMapper);
    selectedActor->GetProperty()->EdgeVisibilityOn();
    selectedActor->GetProperty()->SetColor(
            colors->GetColor3d("Red").GetData());

    selectedActor->GetProperty()->SetLineWidth(5);
    selectedActor->GetProperty()->SetPointSize(5);

    m_vtk_renderer->AddActor(selectedActor);
    m_vtk_window->Render();
}

/*----------------------------------------------------------------------------*/
void GraphicView::update() {
    //means the backend changed, so we update our local pipeline to produce the
    //expected actors

    //clean current actors
    for(auto a: m_vtk_actors) {
        m_vtk_renderer->RemoveActor(a);
        a->Delete();
    }
    if(m_vtk_actors.size() > 2) {
        m_vtk_actors.erase(m_vtk_actors.begin());
        std::cout << "test1" << std::endl;
        m_vtk_actors.erase(m_vtk_actors.begin());
        std::cout << "test2" << std::endl;

    }

    m_vtk_cut_actor->VisibilityOff();
    m_vtk_cut_actor->PickableOff();

    std::cout<<"nb actors left "<<m_vtk_actors.size()<<std::endl;

    m_vtk_actors.clear();
    std::vector<MedusaGrid*> grids = MedusaBackend::getInstance()->grids();
    m_vtk_actors.resize(grids.size());
    int i =0;
    for(auto g:grids) {
        std::cout<<"View: grid name = "<<g->name<<std::endl;
        vtkDataSetMapper *mapper = vtkDataSetMapper::New();
        mapper->SetInputData(g->grid());
        vtkSmartPointer<vtkActor> actor = vtkActor::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetDiffuseColor(0.5,0.5,0.5);
        actor->VisibilityOff();
        m_vtk_actors[i++]=actor;
        m_vtk_renderer->AddActor(actor);
    }

    m_vtk_actors[0]->VisibilityOn();
    m_vtk_actors[0]->GetProperty()->EdgeVisibilityOn();
    m_vtk_actors[0]->GetProperty()->SetLineWidth(0.15);

    if(test == 1){
        m_vtk_renderer->RemoveActor(m_vtk_cut_actor);
        //m_vtk_cut_actor->Delete();
    }
    test = 1;

    if(m_sinGraph_visibility == -1) {
    //    buildSingGraph(MedusaBackend::getInstance()->getSingGraph());
    }
    textMode(0);
    m_vtk_renderer->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::updateMesh(){
    std::vector<MedusaGrid*> grids = MedusaBackend::getInstance()->grids();
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputData(grids[0]->grid());
    m_vtk_actors[0]->SetMapper(mapper);
}
/*----------------------------------------------------------------------------*/
void GraphicView::start() {m_vtk_interactor->Start();}
/*----------------------------------------------------------------------------*/
void GraphicView::createSurface(const std::vector<vtkIdType> AIDs, int ASheetID){

    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();


    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    for(auto id:AIDs)
        ids->InsertNextValue(id);


    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);


    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);


    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());

    std::vector<double> color = {(double)((rand()%256))/256,(double)((rand()%256))/256,(double)((rand()%256))/256};
    m_surf_to_colors.emplace(ASheetID, color);

    vtkDataSetMapper *mapper_selection = vtkDataSetMapper::New();
    mapper_selection->SetInputData(selected);
    mapper_selection->ScalarVisibilityOff();
    vtkSmartPointer<vtkActor> actor_selection = vtkActor::New();
    actor_selection->SetMapper(mapper_selection);
    actor_selection->GetProperty()->EdgeVisibilityOff();
    actor_selection->GetProperty()->SetLineWidth(0.25);
    actor_selection->GetProperty()->SetColor(color[0],color[1],color[2]);

    if(surface_mode == 1){
        actor_selection->VisibilityOff();
    }
    //setOpacityOFF();
    m_vtk_tets_actors.insert(std::pair<vtkSmartPointer<vtkActor>,int>(actor_selection,ASheetID));

    actor_selection->GetProperty()->SetLineWidth(5);
    actor_selection->GetProperty()->SetPointSize(5);
    actor_selection->PickableOff();

    m_vtk_renderer->AddActor(actor_selection);

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::updateDual(int ANbZones) {

    m_nb_zones = ANbZones;


    vtkDataArray *array = MedusaBackend::getInstance()->grids()[0]->grid()->GetCellData()->GetArray("blocks");
    MedusaBackend::getInstance()->grids()[0]->grid()->GetCellData()->SetScalars(array);

    vtkSmartPointer<vtkThreshold> threshold_zone = vtkThreshold::New();
    threshold_zone->SetInputData(MedusaBackend::getInstance()->grids()[0]->grid());
    if(ANbZones == -2){
        threshold_zone->ThresholdByLower(-2);

    }else{
        threshold_zone->ThresholdByUpper(1);
    }
    threshold_zone->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "blocks");
    threshold_zone->Update();
    vtkUnstructuredGrid* thresholdedPolydata = threshold_zone->GetOutput();


    m_mapper_dual->SetInputData(thresholdedPolydata);
    m_mapper_dual->SetScalarModeToUseCellData();
    m_mapper_dual->SetColorModeToMapScalars();

    m_mapper_dual->ScalarVisibilityOn();
    m_mapper_dual->SetScalarRange(-2,10);
    m_mapper_dual->Update();


    m_vtk_actors[0]->SetMapper(m_mapper_dual);

    for(auto a : m_vtk_tets_actors){
        a.first->VisibilityOff();
    }

    m_vtk_renderer->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::updateCut() {

    vtkDataArray *array = MedusaBackend::getInstance()->grids()[0]->grid()->GetCellData()->GetArray("cut");

    MedusaBackend::getInstance()->grids()[0]->grid()->GetCellData()->SetScalars(array);


    vtkSmartPointer<vtkThreshold> threshold_zone = vtkThreshold::New();
    threshold_zone->AddInputData(MedusaBackend::getInstance()->grids()[0]->grid());
    threshold_zone->ThresholdByLower(1);
    threshold_zone->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "cut");
    threshold_zone->Update();
    vtkUnstructuredGrid* thresholdedPolydata = threshold_zone->GetOutput();

    vtkDataSetMapper *m_mapper_cut = vtkDataSetMapper::New();
    m_mapper_cut->SetInputData(thresholdedPolydata);
    m_mapper_cut->SetScalarModeToUseCellData();
    m_mapper_cut->SetColorModeToDefault();
    m_mapper_cut->ScalarVisibilityOff();
    m_mapper_cut->SetScalarRange(0,1);


    m_vtk_actors[0]->SetMapper(m_mapper_cut);
    //m_vtk_actors[0]->GetProperty()->SetColor(0,0,0);

    removeSelection();

    //En fait c'est l'acteur pour les parties coupés
    vtkSmartPointer<vtkThreshold> threshold_cut = vtkThreshold::New();
    threshold_cut->AddInputData(MedusaBackend::getInstance()->grids()[0]->grid());
    threshold_cut->ThresholdByUpper(2);
    threshold_cut->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "cut");
    threshold_cut->Update();
    vtkUnstructuredGrid* thresholdedPolydata_cut = threshold_cut->GetOutput();

    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputData(thresholdedPolydata_cut);
    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarRange(1,2);

    m_vtk_cut_actor->SetMapper(mapper);
    //Attention on a changé l'opacité ici
    m_vtk_cut_actor->VisibilityOff();
    m_vtk_cut_actor->GetProperty()->SetOpacity(0.5);

    m_vtk_renderer->AddActor(m_vtk_cut_actor);



    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::setDualVisibilityON() {
    for(auto a:m_vtk_tets_actors){
        a.first->VisibilityOn();
        a.first->GetProperty()->SetDiffuseColor(0.5,0.5,0.5);
        //a->GetProperty()->SetColor(1.0,0.0,0.0);
    }
    for(auto a:m_vtk_surface_actors){
        a.first->VisibilityOn();
        a.first->GetProperty()->SetDiffuseColor(0.5,0.5,0.5);
        //a->GetProperty()->SetColor(1.0,0.0,0.0);
    }
    MedusaGrid* grid = MedusaBackend::getInstance()->grids()[0];
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputData(grid->grid());
    mapper->SelectColorArray("blocks");
    mapper->SetColorModeToMapScalars();
    mapper->SetScalarModeToUseCellFieldData();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarRange(-1,m_nb_zones-1);
    mapper->Update();


    vtkSmartPointer<vtkActor> actor = vtkActor::New();
    m_vtk_actors[0]->SetMapper(mapper);
    m_vtk_actors[0]->GetProperty()->SetOpacity(1.0);

}
/*----------------------------------------------------------------------------*/
void GraphicView::setDualVisibilityOFF() {
    for(auto a:m_vtk_tets_actors){
        a.first->VisibilityOff();
        //a->GetProperty()->SetColor(1.0,0.0,0.0);
    }
    for(auto a:m_vtk_surface_actors){
        a.first->VisibilityOff();
        //a->GetProperty()->SetColor(1.0,0.0,0.0);
    }

    m_vtk_actors[0]->SetMapper(m_mapper_dual);

    m_vtk_actors[0]->VisibilityOn();
    m_vtk_actors[0]->GetProperty()->SetRepresentationToSurface();
}
/*----------------------------------------------------------------------------*/
void GraphicView::wireframe() {

    if (m_mode == 0) {
        for (auto a:m_vtk_actors) {
            a->GetProperty()->SetRepresentationToSurface();
        }
        m_vtk_actors[0]->GetProperty()->SetRepresentationToWireframe();
    }else{
        for (auto a:m_vtk_actors) {
            a->GetProperty()->SetRepresentationToWireframe();
        }
        m_vtk_actors[0]->GetProperty()->SetRepresentationToSurface();
    }
    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::wireframeAll(){

    for(auto a:m_vtk_actors){
        a->GetProperty()->SetRepresentationToWireframe();
    }
    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::solid(){
    for (auto a:m_vtk_actors) {
        a->GetProperty()->SetRepresentationToSurface();
    }
    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::opacity(){

    if(m_opacity){
        setOpacityOFF();
    }else{
        setOpacityON();
    }

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::setOpacityON(){
    if(m_mode == 0){
        m_vtk_actors[0]->GetProperty()->SetOpacity(1.0);
    }else{
        setDualVisibilityOFF();
    }
    m_opacity = true;
}
/*----------------------------------------------------------------------------*/
void GraphicView::setOpacityOFF(){
    if(m_mode == 0){
        m_vtk_actors[0]->GetProperty()->SetOpacity(0.25);
    }else{
        setDualVisibilityON();
    }
    m_opacity = false;
}
/*----------------------------------------------------------------------------*/
void GraphicView::toggleCutOpacity(){
    if(m_op_cut == 2){
        m_vtk_cut_actor->GetProperty()->SetOpacity(0);
        m_op_cut = 0;
    }else if(m_op_cut == 0){
        m_vtk_cut_actor->GetProperty()->SetOpacity(0.5);
        m_op_cut = 1;
    }else{
        m_vtk_cut_actor->GetProperty()->SetOpacity(1);
        m_op_cut = 2;
    }
}
/*----------------------------------------------------------------------------*/
// Maybe not necessary
void GraphicView::unlockSheet(){
    if(surface_mode == 0) {
        for (auto a:m_vtk_tets_actors) {
            a.first->PickableOn();
        }
    }else{
        for (auto a:m_vtk_surface_actors) {
            a.first->PickableOn();
        }
    }
    m_vtk_actors[0]->PickableOff();
}
/*----------------------------------------------------------------------------*/
// Maybe not necessary
void GraphicView::unlockLines(){

    for (auto a:m_vtk_sing_actors) {
        std::cout<<a.first<<std::endl;
        a.first->PickableOn();
    }

    m_vtk_actors[0]->PickableOff();
}
/*----------------------------------------------------------------------------*/
int GraphicView::selectLine(const vtkSmartPointer<vtkActor> ADeletedActor){

    ADeletedActor->GetProperty()->SetDiffuseColor(0,0,1);
    for(auto const &line : m_vtk_sing_actors){
        if(ADeletedActor == line.first){
            return line.second;
        }
    }
}
/*----------------------------------------------------------------------------*/
void GraphicView::showSetOfTets(std::vector<vtkIdType> AIDs){

    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    for(auto id:AIDs)
        ids->InsertNextValue(id);


    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);


    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);


    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());


    vtkSmartPointer<vtkDataSetMapper> mapper = vtkDataSetMapper::New();
    mapper->SetInputData(selected);
    mapper->ScalarVisibilityOff();

    m_vtk_tets_set->SetMapper(mapper);
    m_vtk_tets_set->VisibilityOn();
    m_vtk_tets_set->GetProperty()->EdgeVisibilityOff();
    m_vtk_tets_set->GetProperty()->SetColor(1,0,0);
    m_vtk_tets_set->GetProperty()->SetLineWidth(5);
    m_vtk_tets_set->GetProperty()->SetPointSize(5);
    m_vtk_tets_set->PickableOff();


    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::addTetSetActor(){
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(m_vtk_tets_set->GetMapper());
    m_vtk_sets_actors.emplace_back(actor);
}
/*----------------------------------------------------------------------------*/
void GraphicView::removeLastSelection(){
    auto actor = m_vtk_selection_actors.back();
    m_vtk_selection_actors.pop_back();
    m_vtk_renderer->RemoveActor(actor);
    actor->Delete();
}
/*----------------------------------------------------------------------------*/
void GraphicView::removeSelection(){
    for(int i = 0; i<m_vtk_selection_actors.size();i++){
        m_vtk_renderer->RemoveActor(m_vtk_selection_actors[i]);
        m_vtk_selection_actors[i]->Delete();
    }
}
/*----------------------------------------------------------------------------*/
int GraphicView::removeActor(const vtkSmartPointer<vtkActor> ADeletedActor){

    int sheet_id = 0;
    if(surface_mode == 0) {
        for (auto a:m_vtk_tets_actors) {
            a.first->PickableOff();
            if (a.first == ADeletedActor) {
                sheet_id = a.second;
                m_vtk_renderer->RemoveActor(a.first);
                a.first->Delete();
            }
        }
        for(auto a:m_vtk_surface_actors){
            if(a.second == sheet_id){
                m_vtk_renderer->RemoveActor(a.first);
                a.first->Delete();
            }
        }
    } else{
        for (auto a:m_vtk_surface_actors) {
            a.first->PickableOff();
            if (a.first == ADeletedActor) {
                sheet_id = a.second;
                m_vtk_renderer->RemoveActor(a.first);
                a.first->Delete();
            }
        }
        for(auto a:m_vtk_tets_actors){
            if(a.second == sheet_id){
                m_vtk_renderer->RemoveActor(a.first);
                a.first->Delete();
            }
        }
    }
    m_surf_to_colors.erase(sheet_id);
    m_vtk_window->Render();
    m_vtk_actors[0]->PickableOn();
    return sheet_id;
}
/*----------------------------------------------------------------------------*/
void GraphicView::removeActor(int ADeletedID){

    for(auto a:m_vtk_tets_actors){
        if(a.second == ADeletedID){
            vtkSmartPointer<vtkActor> actor = a.first;
            m_vtk_renderer->RemoveActor(actor);
            a.first->Delete();
        }
    }
    for(auto a:m_vtk_surface_actors){
        if(a.second == ADeletedID){
            vtkSmartPointer<vtkActor> actor = a.first;
            m_vtk_renderer->RemoveActor(actor);
            a.first->Delete();
        }
    }
    m_surf_to_colors.erase(ADeletedID);
    m_vtk_window->Render();
    m_vtk_actors[0]->PickableOn();
}
/*----------------------------------------------------------------------------*/
void GraphicView::showAxis(double ACoords[3][3]){

// Create the polydata where we will store all the geometric data
    vtkSmartPointer<vtkPolyData> linesPolyData =
            vtkSmartPointer<vtkPolyData>::New();

    double resize = (std::sqrt((m_vtk_window->GetSize()[0])*(m_vtk_window->GetSize()[0])+(m_vtk_window->GetSize()[1])*(m_vtk_window->GetSize()[1])))*0.01;

    // Create a vtkPoints container and store the points in it
    vtkSmartPointer<vtkPoints> pts =
            vtkSmartPointer<vtkPoints>::New();
    pts->InsertNextPoint(ACoords[0]);
    pts->InsertNextPoint(ACoords[1]);
    pts->InsertNextPoint(ACoords[2]);


    // Add the points to the polydata container
    linesPolyData->SetPoints(pts);


    // Create the first line (between Origin and P0)
    vtkSmartPointer<vtkLine> line0 =
            vtkSmartPointer<vtkLine>::New();
    line0->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
    line0->GetPointIds()->SetId(1, 1); // the second 1 is the index of P0 in linesPolyData's points

    // Create the second line (between Origin and P1)
    vtkSmartPointer<vtkLine> line1 =
            vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
    line1->GetPointIds()->SetId(1, 2); // 2 is the index of P1 in linesPolyData's points

    // Create a vtkCellArray container and store the lines in it
    vtkSmartPointer<vtkCellArray> lines =
            vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(line0);
    lines->InsertNextCell(line1);

    // Add the lines to the polydata container
    linesPolyData->SetLines(lines);


    // Create two colors - one for each line
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};

    // Create a vtkUnsignedCharArray container and store the colors in it
    vtkSmartPointer<vtkUnsignedCharArray> colors =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->InsertNextTypedTuple(red);
    colors->InsertNextTypedTuple(green);

    // Color the lines.
    // SetScalars() automatically associates the values in the data array passed as parameter
    // to the elements in the same indices of the cell data array on which it is called.
    // This means the first component (red) of the colors array
    // is matched with the first component of the cell array (line 0)
    // and the second component (green) of the colors array
    // is matched with the second component of the cell array (line 1)
    linesPolyData->GetCellData()->SetScalars(colors);


    // Setup the visualization pipeline
    vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(linesPolyData);

    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(4);

    m_vtk_renderer->AddActor(actor);
    m_vtk_window->Render();

}
/*----------------------------------------------------------------------------*/
void GraphicView::removeAxis(){
    vtkSmartPointer<vtkActor> last_act = m_vtk_renderer->GetActors()->GetLastItem();
    m_vtk_renderer->RemoveActor(last_act);
    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::setModeDomain(){
    m_mode = 0;
    setDualVisibilityON();
    m_vtk_actors[0]->GetMapper()->ScalarVisibilityOff();
    m_vtk_actors[0]->GetProperty()->SetDiffuseColor(0.5,0.5,0.5);
    for(auto a : m_vtk_tets_actors){
        std::vector<double> colors;
        colors.resize(3);
        colors = m_surf_to_colors[a.second];
        a.first->GetProperty()->SetDiffuseColor(colors[0],colors[1],colors[2]);
        a.first->VisibilityOn();
    }
    for(auto a : m_vtk_surface_actors){
        std::vector<double> colors;
        colors.resize(3);
        colors = m_surf_to_colors[a.second];
        a.first->GetProperty()->SetDiffuseColor(colors[0],colors[1],colors[2]);
        a.first->VisibilityOn();
    }
    refresh();
}
/*----------------------------------------------------------------------------*/
void GraphicView::setModeDual(){
    m_mode = 1;
}
/*----------------------------------------------------------------------------*/
void GraphicView::viewBlocks() {
    //update();
    m_vtk_actors[0]->VisibilityOff();
    m_vtk_actors[0]->GetMapper()->SetScalarModeToUseCellData();
    m_vtk_actors[0]->GetMapper()->SetColorModeToMapScalars();
    m_vtk_actors[0]->GetMapper()->ScalarVisibilityOn();

    for(auto a : m_vtk_actors){
        a->VisibilityOff();
    }

    m_vtk_actors[1]->VisibilityOn();
    m_vtk_actors[1]->GetProperty()->SetRepresentationToSurface();
    m_vtk_actors[1]->GetProperty()->EdgeVisibilityOn();

    m_block_mode = true;
}
/*----------------------------------------------------------------------------*/
void GraphicView::viewSurface() {
    if(surface_mode == 0) {
        for(auto a : m_vtk_tets_actors){
            a.first->VisibilityOff();
        }
        for(auto a : m_vtk_surface_actors){
            a.first->VisibilityOn();
        }
        surface_mode = 1;
    }else{
        for(auto a : m_vtk_tets_actors){
            a.first->VisibilityOn();
        }
        for(auto a : m_vtk_surface_actors){
            a.first->VisibilityOff();
        }
        surface_mode = 0;
    }

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::createPlanSurface(int ASheetID, MedusaGrid* ASurf) {

    double colors[3];

    for(auto c : m_surf_to_colors){
        if(c.first == ASheetID){
            colors[0] = c.second[0];
            colors[1] = c.second[1];
            colors[2] = c.second[2];
        }
    }

    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputData(ASurf->grid());
    vtkSmartPointer<vtkActor> actor = vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetDiffuseColor(colors);
    m_vtk_surface_actors.insert(std::pair<vtkSmartPointer<vtkActor>,int>(actor,ASheetID));
    if(surface_mode == 0){
        actor->VisibilityOff();
    }

    m_vtk_renderer->AddActor(actor);

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::toggleBlocks() {
    if(!m_block_mode){
        m_vtk_actors[0]->VisibilityOff();
        m_vtk_actors[1]->VisibilityOn();
        m_block_mode = true;
    } else{
        m_vtk_actors[0]->VisibilityOn();
        m_vtk_actors[1]->VisibilityOff();
        m_block_mode = false;
    }
}
/*----------------------------------------------------------------------------*/
void GraphicView::buildSingGraph(std::vector<vtkIdType> AIDs){

    if(AIDs.empty()) {
        return;
    }

    m_sinGraph_visibility = 1;

    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();


    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    for(auto id:AIDs)
        ids->InsertNextValue(id);


    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);


    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);


    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());



    vtkDataSetMapper *mapper_selection = vtkDataSetMapper::New();
    mapper_selection->SetInputData(selected);
    vtkSmartPointer<vtkActor> actor_selection = vtkActor::New();
    actor_selection->SetMapper(mapper_selection);
    actor_selection->GetProperty()->EdgeVisibilityOff();
    actor_selection->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

    actor_selection->GetProperty()->SetLineWidth(5);
    actor_selection->GetProperty()->SetPointSize(5);
    actor_selection->PickableOff();

    m_singGraph = actor_selection;

    m_vtk_actors.push_back(actor_selection);
    m_vtk_renderer->AddActor(actor_selection);

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::resetSingLine(){
    for(auto line : m_vtk_sing_actors){
        m_vtk_renderer->RemoveActor(line.first);
        line.first->Delete();
    }
    m_vtk_sing_actors.clear();
}
/*----------------------------------------------------------------------------*/
void GraphicView::buildSingLine(int AID, std::vector<vtkIdType> AIDs) {

    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();

    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    for(auto id:AIDs)
        ids->InsertNextValue(id);

    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);


    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);


    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());



    vtkDataSetMapper *mapper_selection = vtkDataSetMapper::New();
    mapper_selection->SetInputData(selected);
    vtkSmartPointer<vtkActor> actor_selection = vtkActor::New();
    actor_selection->SetMapper(mapper_selection);
    actor_selection->GetProperty()->EdgeVisibilityOff();
    actor_selection->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

    actor_selection->GetProperty()->SetLineWidth(5);
    actor_selection->GetProperty()->SetPointSize(5);
    actor_selection->PickableOff();



    m_vtk_sing_actors.emplace(actor_selection, AID);
    m_vtk_renderer->AddActor(actor_selection);

    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::toggleSingGraph(){

    if(m_sinGraph_visibility){
        m_singGraph->VisibilityOff();
        m_sinGraph_visibility = false;
    }else{
        m_singGraph->VisibilityOn();
        m_sinGraph_visibility = true;
    }
}
/*----------------------------------------------------------------------------*/
void GraphicView::surfacePicked(std::vector<vtkIdType> AIDs) {
    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();


    vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    for(auto id:AIDs)
        ids->InsertNextValue(id);


    vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::CELL);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);


    vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);


    vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, MedusaBackend::getInstance()->grids()[0]->grid());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // In selection
    vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());



    vtkDataSetMapper *mapper_selection = vtkDataSetMapper::New();
    mapper_selection->SetInputData(selected);
    vtkSmartPointer<vtkActor> actor_selection = vtkActor::New();
    actor_selection->SetMapper(mapper_selection);
    actor_selection->GetProperty()->EdgeVisibilityOff();
    actor_selection->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

    actor_selection->GetProperty()->SetLineWidth(5);
    actor_selection->GetProperty()->SetPointSize(5);
    actor_selection->PickableOff();

    m_vtk_selection_actors.push_back(actor_selection);

    m_vtk_renderer->AddActor(actor_selection);

    m_vtk_window->Render();

}
/*----------------------------------------------------------------------------*/
void GraphicView::toggleBRep(){
    //m_vtk_actors[0]->VisibilityOff();
    m_vtk_actors[1]->VisibilityOn();
    m_vtk_actors[1]->PickableOff();
    //m_vtk_actors[0]->PickableOn();

    MedusaGrid* g = MedusaBackend::getInstance()->grids()[1];


    vtkDataArray *array = g->grid()->GetCellData()->GetArray("BND_Color");

    MedusaBackend::getInstance()->grids()[1]->grid()->GetCellData()->SetScalars(array);

    int nbSurfaces = MedusaBackend::getInstance()->grids()[1]->getNbGeomSurfs();

    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputData(g->grid());
    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarRange(0,nbSurfaces);

    mapper->Update();

    m_vtk_actors[1]->SetMapper(mapper);

}
/*----------------------------------------------------------------------------*/
void GraphicView::textMode(int AMode){
    // Setup the text and add it to the renderer

    switch (AMode){
        case 0:
            m_textActor->SetInput ( "Input mode" );
            break;
        case 1:
            m_textActor->SetInput ( "Frame mode" );
            break;
        case 2:
            m_textActor->SetInput ( "Geom mode" );
            break;
        case 3:
            m_textActor->SetInput ( "Dublo mode" );
            break;
    }
    m_vtk_window->Render();
}
/*----------------------------------------------------------------------------*/
void GraphicView::hideSetOfTets() {
    m_vtk_tets_set->VisibilityOff();
}
/*----------------------------------------------------------------------------*/
void GraphicView::showPoints(std::vector<std::vector<double>> coords){
    vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();

    for(auto coord : coords){
        vtkSmartPointer<vtkSphereSource> sphereSource =
                vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(coord[0],coord[1],coord[2]);
        sphereSource->SetRadius(0.25);


        vtkSmartPointer<vtkPolyDataMapper> mapper =
                vtkSmartPointer<vtkPolyDataMapper>::New();

        mapper->SetInputConnection(sphereSource->GetOutputPort());
        vtkSmartPointer<vtkActor> actor =
                vtkActor::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

        m_vtk_renderer->AddActor(actor);
        m_vtk_actors.push_back(actor);
    }

    m_vtk_window->Render();



}
/*----------------------------------------------------------------------------*/
void GraphicView::resetDual() {
    for(auto a : m_vtk_tets_actors){
        m_vtk_renderer->RemoveActor(a.first);
        a.first->Delete();
    }

    m_vtk_tets_actors.clear();

    for(auto a : m_vtk_surface_actors){
        m_vtk_renderer->RemoveActor(a.first);
        a.first->Delete();
    }

    m_vtk_surface_actors.clear();

}
