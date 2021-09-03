/*----------------------------------------------------------------------------*/
#include <medusa/Medusa.h>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout<<"Medusa (Mesh Blocking)"<<std::endl;

    //==================================================================
    std::string fIn;

    if (argc < 2)
        //throw gmds::GMDSException("No file loaded");
        std::cout<<"No file loaded "<<std::endl;

    if (argc >= 2) {
        fIn = std::string(argv[1]);
        std::cout << "INPUT FILE: " << fIn << std::endl;
    }

    Medusa m(Medusa::DUAL_BLOCKING);
    m.newGraphicView("Medusa 1", GraphicView::View3D);
    m.newTextView();
    if(fIn.length() == 0){
        std::cout<<"INPUT FILE : "<<std::endl;
        std::cin >> fIn;
    }
    m.setMode(Medusa::GEOMETRY);
    m.load(fIn);
    m.launch();



    return 0;
}
/*----------------------------------------------------------------------------*/

