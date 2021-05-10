/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDUSA_H
#define GMDS_MEDUSA_H
/*----------------------------------------------------------------------------*/
#include <medusa/view/GraphicView.h>
#include <medusa/view/TextView.h>
#include <medusa/control/MediatorControl.h>
/*----------------------------------------------------------------------------*/
namespace medusa {
    class Medusa {
    public:
        enum Mode{DUAL_BLOCKING, EXPLORATION};
        Medusa(const Mode AMode);
        void setMode(const Mode AMode);

        virtual ~Medusa();
        void load(const std::string& AFileName);
        void launch();
        GraphicView* newGraphicView(const std::string& AName, GraphicView::ViewType AType);
        TextView* newTextView();
    private:
        void init(const Medusa::Mode AMode);
    private:
        GraphicView* m_last_gview;
        std::vector<GraphicView*> m_gviews;
        std::vector<TextView*> m_tviews;


    };
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_MEDUSA_H
/*----------------------------------------------------------------------------*/
