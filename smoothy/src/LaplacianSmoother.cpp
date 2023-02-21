/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/LaplacianSmoother.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
LaplacianSmoother::LaplacianSmoother(GeomMeshLinker* ALinker)
: AbstractSmoother(ALinker) {}
/*----------------------------------------------------------------------------*/
void LaplacianSmoother::smoothCurves(const int ANbIterations) {

    std::vector<GeomCurve*> curves;
    m_linker->geometry()->getCurves(curves);
    Mesh* m = m_linker->mesh();

    for(auto i=0;i<ANbIterations;i++) {
        for (auto c:curves) {

            //We get all the nodes of the curve
            auto ids_of_nodes_on_curve = m_c2n[c->id()];
            //We smooth each node in a naive way
            for (auto nid: ids_of_nodes_on_curve) {
                Node ni = m->get<Node>(nid);
                math::Point pi = ni.point();
                std::vector<TCellID> adj_n = m_n2n[nid];
                math::Point pnew(0, 0, 0);
                if(adj_n.size()>2) {
                    throw GMDSException("LaplacianSmoother error: a node on curve has more that 2 neighhbors.");
                }
                for (auto adj_id : adj_n) {
                    math::Point pj = m->get<Node>(adj_id).point();
                    pnew = pnew+ pj;
                }
                if(adj_n.size()!=0){
                    pnew =(1.0 / adj_n.size())*pnew;
                    pi = 0.75*pi+0.25*pnew;
                    c->project(pi);
                    ni.setPoint(pi);
                }
            }
        }
    }

}
/*----------------------------------------------------------------------------*/
void LaplacianSmoother::smoothSurfaces(const int ANbIterations) {


    std::vector<GeomSurface*> surfs;
    m_linker->geometry()->getSurfaces(surfs);
    Mesh* m = m_linker->mesh();
    for(auto i=0;i<ANbIterations;i++) {
        for (auto s:surfs) {
            //We get all the nodes of the surface
            auto ids_of_nodes_on_surf = m_s2n[s->id()];
            //We smooth each node in a naive way
            for (auto nid: ids_of_nodes_on_surf) {
                Node ni = m->get<Node>(nid);
                std::vector<TCellID> adj_n = m_n2n[nid];
                math::Point pi(0, 0, 0);
                for (auto adj_id : adj_n) {
                    math::Point pj = m->get<Node>(adj_id).point();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());
                try{
                    s->project(pi);

                }
                catch(GMDSException& e){}
                ni.setPoint(pi);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void LaplacianSmoother::smoothVolumes(const int ANbIterations) {

    Mesh* m = m_linker->mesh();

    for(auto i=0;i<ANbIterations;i++) {
        for (auto vol_info:m_v2n) {

            //We get all the nodes of the first volume (or void)
            auto ids_of_nodes_in_vol = vol_info.second;
            //We smooth each node in a naive way
            for (auto nid: ids_of_nodes_in_vol) {
                Node ni = m->get<Node>(nid);
                std::vector<TCellID> adj_n = m_n2n[nid];
                math::Point pi(0, 0, 0);
                for (auto adj_id : adj_n) {
                    math::Point pj = m->get<Node>(adj_id).point();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());
                ni.setPoint(pi);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
