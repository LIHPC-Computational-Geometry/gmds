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
        throw gmds::GMDSException("Wrong parameters");

    if (argc >= 2) {
        fIn = std::string(argv[1]);
        std::cout << "INPUT FILE: " << fIn << std::endl;
    }

    Medusa m(Medusa::DUAL_BLOCKING);
    m.newGraphicView("Medusa 1", GraphicView::View3D);
    m.newTextView();
    m.load(fIn);
    m.launch();


    return 0;
}
/*----------------------------------------------------------------------------*/

