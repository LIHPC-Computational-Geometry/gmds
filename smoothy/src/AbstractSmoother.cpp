/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/AbstractSmoother.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
AbstractSmoother::AbstractSmoother(GeomMeshLinker* ALinker)
:m_linker(ALinker)
{
    init();
}
/*----------------------------------------------------------------------------*/
bool AbstractSmoother::isValid() const
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
void AbstractSmoother::init() {

    //NO I NEED EDGE ON CURVE AND FACE ON SURFACE. OTHERWISE, I CANNOT DO SOMETHING
    //GOOD FOR THIN AREA, WHICH WILL OCCUR VERY OFTEN IN FACT

    Mesh *m = m_linker->mesh();
    for (auto n_id:m->nodes()) {
        auto geom_info = m_linker->getGeomInfo<Node>(n_id);
        if (geom_info.first == GeomMeshLinker::LinkCurve) {
            //Means on a curve
            m_c2n[geom_info.second].push_back(n_id);
        } else if (geom_info.first == GeomMeshLinker::LinkSurface) {
            //Means on a surface
            m_s2n[geom_info.second].push_back(n_id);
        } else if (geom_info.first == GeomMeshLinker::NoLink ||
        geom_info.first == GeomMeshLinker::LinkVolume) {
            //Means in a volume
            m_v2n[geom_info.second].push_back(n_id);
        }

    }

    //Now we build the local N2N connection that is useful for the smoothing
    for (auto n_id:m->nodes()) {
        auto geom_dim = m_linker->getGeomDim<Node>(n_id);
        auto geom_id = m_linker->getGeomId<Node>(n_id);

        if (geom_dim == GeomMeshLinker::LinkCurve) {
            //Means on a curve
            std::vector<Edge> edges = m->get<Node>(n_id).get<Edge>();

            for (auto e:edges) {
                if (m_linker->getGeomDim(e) == GeomMeshLinker::LinkCurve &&
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
        } else if (geom_dim == GeomMeshLinker::LinkSurface) {
            //Means on a surface
            std::vector<Face> faces = m->get<Node>(n_id).get<Face>();

            for (auto f:faces) {
                if (m_linker->getGeomDim(f) == GeomMeshLinker::LinkSurface &&
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
        } else if (geom_dim == GeomMeshLinker::NoLink ||
        geom_dim == GeomMeshLinker::LinkVolume) {
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
