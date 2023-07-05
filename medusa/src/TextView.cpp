/*----------------------------------------------------------------------------*/
#include <medusa/view/TextView.h>
#include <gmds/utils/Exception.h>
#include <iostream>
#include <medusa/model/MedusaBackEnd.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
TextView::TextView(){}
/*----------------------------------------------------------------------------*/
void TextView::refresh() {
    throw gmds::GMDSException("Not yet implemented");

}/*----------------------------------------------------------------------------*/
void TextView::update() {
    std::cout<<"Update of text view"<<std::endl;

}
/*----------------------------------------------------------------------------*/
void TextView::selectVTKCell(int ADIm, vtkIdType ACellID) {
    std::vector<MedusaGrid*> grids = MedusaBackend::getInstance()->grids();
    gmds::TCellID gmds_id = grids[0]->getGMDSCellID(ADIm,ACellID);
    std::cout<<"Selection ("<<ADIm<<", "<<gmds_id<<")"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void TextView::createSurface(const std::vector<vtkIdType> AIDs, int ASheetID) {
    std::cout<<"Nb tet in dual surface "<<ASheetID<<": "<<AIDs.size()<<std::endl;
}
/*----------------------------------------------------------------------------*/
void TextView::help() {
    std::cout<<"w - Wireframe"<<std::endl;
    std::cout<<"s - Solid"<<std::endl;
    std::cout<<"p - Pick/create surf"<<std::endl;
    std::cout<<"r - Remove"<<std::endl;
    std::cout<<"h - Help"<<std::endl;
    std::cout<<"g - Generate dual"<<std::endl;
    std::cout<<"b - Block generation"<<std::endl;
    std::cout<<"v - toggle surface View"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void TextView::generateDual() {
    std::cout<<"Dual generated"<<std::endl;
}
/*----------------------------------------------------------------------------*/
void TextView::dualError() {
    std::cout<<"ERROR: The dual is not generated yet"<<std::endl;
}
/*----------------------------------------------------------------------------*/
int TextView::getInputAxis() {
    char axis = '0';
    while(true){
        std::cout<<"Select dual surface normal direction : '1' for red, '2' for green or '3' to create a boundary surface"<<std::endl;
        cin >> axis;
        switch(axis){
            case '1':
                return 1;
            case '2':
                return 2;
            case '3':
                return 3;
            case 'q':
                return 0;
            default:
                std::cout<<"Wrong input, please retry."<<std::endl;
                cin.clear();
        }
    }
}
/*----------------------------------------------------------------------------*/
void TextView::remove(int ASheetID) {
    std::cout<<"The dual surface "<<ASheetID<<" has been removed"<<std::endl;
}
/*----------------------------------------------------------------------------*/
