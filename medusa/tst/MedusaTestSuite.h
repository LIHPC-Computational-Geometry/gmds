/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include <medusa/Medusa.h>
/*----------------------------------------------------------------------------*/
using namespace medusa;
/*----------------------------------------------------------------------------*/
TEST(MedusaTest, init) {
    Medusa m(Medusa::EXPLORATION);
    m.newGraphicView("Medusa 1", GraphicView::View3D);
    m.newTextView();
}
