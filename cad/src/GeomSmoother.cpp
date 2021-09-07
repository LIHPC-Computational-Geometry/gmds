/*----------------------------------------------------------------------------*/
#include "gmds/cad/GeomSmoother.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
GeomSmoother::GeomSmoother(GeomMeshLinker* ALinker)
:m_linker(ALinker)
{
    init();
}
/*----------------------------------------------------------------------------*/
bool GeomSmoother::isValid() const
{
    bool valid_input=true;
    MeshModel model = m_linker->mesh()->getModel();
    if (!model.has(R)   ||
    !model.has(F)   ||
    !model.has(E)   ||
    !model.has(N)   ||
    !model.has(R2N) ||
    !model.has(F2N) ||
    !model.has(E2N) ||
    !model.has(N2E) ||
    !model.has(N2F) )
        valid_input= false;

    if(!valid_input)
        return false;


    return true;
}
/*----------------------------------------------------------------------------*/
void GeomSmoother::init() {

    //NO I NEED EDGE ON CURVE AND FACE ON SURFACE. OTHERWISE, I CANNOT DO SOMETHING
    //GOOD FOR THIN AREA, WHICH WILL OCCUR VERY OFTEN IN FACT


    Mesh *m = m_linker->mesh();
    for (auto n_id:m->nodes()) {
        auto geom_info = m_linker->getGeomInfo<Node>(n_id);

        if (geom_info.first == GeomMeshLinker::LINK_CURVE) {
            //Means on a curve
            m_c2n[geom_info.second].push_back(n_id);
        } else if (geom_info.first == GeomMeshLinker::LINK_SURFACE) {
            //Means on a surface
            m_s2n[geom_info.second].push_back(n_id);
        } else if (geom_info.first == GeomMeshLinker::NO_LINK ||
        geom_info.first == GeomMeshLinker::LINK_VOLUME) {
            //Means in a volume
            m_v2n[geom_info.second].push_back(n_id);
        }

    }

    //Now we build the local N2N connection that is useful for the smoothing
    for (auto n_id:m->nodes()) {
        auto geom_dim = m_linker->getGeomDim<Node>(n_id);
        auto geom_id = m_linker->getGeomId<Node>(n_id);

        if (geom_dim == GeomMeshLinker::LINK_CURVE) {

            //Means on a curve
            std::vector<Edge> edges = m->get<Node>(n_id).get<Edge>();

            for (auto e:edges) {
                if (m_linker->getGeomDim(e) == GeomMeshLinker::LINK_CURVE &&
                m_linker->getGeomId(e) == geom_id) {
                    //Edge e is linked to the same curve as n_id
                    std::vector<TCellID> nids = e.getIDs<Node>();
                    if (nids[0] == n_id) {
                        m_n2n[n_id].push_back(nids[1]);
                    } else {
                        m_n2n[n_id].push_back(nids[0]);
                    }
                }
            }
        } else if (geom_dim == GeomMeshLinker::LINK_SURFACE) {
            //Means on a surface
            std::vector<Face> faces = m->get<Node>(n_id).get<Face>();

            for (auto f:faces) {
                if (m_linker->getGeomDim(f) == GeomMeshLinker::LINK_SURFACE &&
                m_linker->getGeomId(f) == geom_id) {
                    //Face f is linked to the same surface as n_id
                    std::vector<TCellID> nids = f.getIDs<Node>();
                    for (auto ni:nids) {
                        if (ni != n_id) {
                            m_n2n[n_id].push_back(ni);
                        }
                    }
                }
            }
        } else if (geom_dim == GeomMeshLinker::NO_LINK ||
        geom_dim == GeomMeshLinker::LINK_VOLUME) {
            //Means in a volume
            std::vector<Region> regions = m->get<Node>(n_id).get<Region>();

            for (auto r:regions) {
                std::vector<TCellID> nids = r.getIDs<Node>();
                for (auto ni:nids) {
                    if (ni != n_id) {
                        m_n2n[n_id].push_back(ni);
                    }
                }
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void GeomSmoother::smoothCurves(const int ANbIterations) {

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
                std::vector<TCellID> adj_n = m_n2n[nid];
                math::Point pi(0, 0, 0);
                for (auto adj_id : adj_n) {
                    math::Point pj = m->get<Node>(adj_id).getPoint();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());

                c->project(pi);
                ni.setPoint(pi);
            }
        }
    }

}
/*----------------------------------------------------------------------------*/
void GeomSmoother::smoothSurfaces(const int ANbIterations) {


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
                    math::Point pj = m->get<Node>(adj_id).getPoint();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());
                s->project(pi);
                ni.setPoint(pi);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void GeomSmoother::smoothVolumes(const int ANbIterations) {

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
                    math::Point pj = m->get<Node>(adj_id).getPoint();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());
                ni.setPoint(pi);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
