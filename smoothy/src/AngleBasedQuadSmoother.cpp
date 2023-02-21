/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/AngleBasedQuadSmoother.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
AngleBasedQuadSmoother::AngleBasedQuadSmoother(GeomMeshLinker* ALinker)
:AbstractSmoother(ALinker)
{}
/*----------------------------------------------------------------------------*/
void AngleBasedQuadSmoother::smooth(const int ANbIterations) {
    //we iteratively smooth curves and surfaces
    for (auto i = 0; i < ANbIterations; i++) {
        smoothSurfaces();
        smoothCurves();
    }
}
/*----------------------------------------------------------------------------*/
bool AngleBasedQuadSmoother::isValid() const {
    bool model_valid = AbstractSmoother::isValid();
    if(!model_valid)
        return false;
    for(auto f_id:m_linker->mesh()->faces()){
        if(m_linker->mesh()->get<Face>(f_id).type()!=GMDS_QUAD)
            return false;
    }
    return true;
}

/*----------------------------------------------------------------------------*/
void AngleBasedQuadSmoother::smoothCurves() {
    std::vector<GeomCurve*> curves;
    m_linker->geometry()->getCurves(curves);
    Mesh* m = m_linker->mesh();

    for (auto c:curves) {
        //We get all the nodes of the curve
        auto ids_of_nodes_on_curve = m_c2n[c->id()];

        //we get adjacent surface ids
        std::vector<GeomSurface*> adj_surf = c->surfaces();
        //We smooth each node
        for (auto nid: ids_of_nodes_on_curve) {
            Node ni = m->get<Node>(nid);
            math::Point pi = ni.point();
            std::vector<TCellID> adj_n = m_n2n[nid];
            if(adj_n.size()!=2) {
                std::cout<<"Smoother error. Node "<<ni.id()<<" located at "<<ni.point()<<" has the following neighbours ";
                for(auto a : adj_n){
                    std::cout<< a<<" ";
                }
                std::cout<<std::endl;
                throw GMDSException("Smoother error: a node on curve has more that 2 neighhbors.");
            }
            //we compute the location we would like for the curve only
            math::Point curve_contribution(0, 0, 0);
            for (auto adj_id : adj_n) {
                math::Point pj = m->get<Node>(adj_id).point();
                curve_contribution = curve_contribution + pj;
            }
            curve_contribution = (1.0 / adj_n.size()) * curve_contribution;

            math::Segment curve_linear_approx(m->get<Node>(adj_n[0]).point(),
                                              m->get<Node>(adj_n[1]).point());

            //For each surface adjacent to the curve, we check if we have to try to compute orthogonal edges
            std::vector<math::Point> surfaces_contribution;
            std::vector<TCellID> adj_face_ids = ni.getIDs<Face>();
            for(auto s:adj_surf){
                auto s_id = s->id();
                std::vector<TCellID> faces_on_surf;
                for(auto f_id:adj_face_ids) {

                    if (m_linker->getGeomDim<Face>(f_id) == GeomMeshLinker::LinkSurface &&
                        m_linker->getGeomId<Face>(f_id) ==s_id){
                      faces_on_surf.push_back(f_id);
                  }
                }
                if(faces_on_surf.size()==2){
                    //means those faces share an edge that may be orthogonal to the curve
                    Face f0 = m_linker->mesh()->get<Face>(faces_on_surf[0]);
                    Face f1 = m_linker->mesh()->get<Face>(faces_on_surf[1]);
                    std::vector<TCellID > nodes_f01 = m_linker->mesh()->getCommonNodes(f0,f1);
                    TCellID opp_node_id = (nodes_f01[0]==ni.id())?nodes_f01[1]:nodes_f01[0];
                    surfaces_contribution.push_back(curve_linear_approx.project(
                            m_linker->mesh()->get<Node>(opp_node_id).point()));
                }
            }
            math::Point geom_contribution =curve_contribution;
            for(auto surf_point:surfaces_contribution){
                geom_contribution = geom_contribution+surf_point;
            }
            geom_contribution = (1.0/(1.0+surfaces_contribution.size()))*geom_contribution;
            pi = geom_contribution;
            c->project(pi);
            ni.setPoint(pi);
        }
    }

}
/*----------------------------------------------------------------------------*/
void AngleBasedQuadSmoother::smoothSurfaces() {

    std::vector<GeomSurface*> surfs;
    m_linker->geometry()->getSurfaces(surfs);
    Mesh* m = m_linker->mesh();
    for (auto s:surfs) {
        //We get all the nodes of the surface
        auto ids_of_nodes_on_surf = m_s2n[s->id()];
        //We smooth each node in a laplacian way when its number of adjacent
        // faces is not equal to 4.
        // When it is 4? we used the equipotential smoothing,
        // i.e. scheme C of https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.79.9662&rep=rep1&type=pdf

        for (auto nid: ids_of_nodes_on_surf) {
            Node ni = m->get<Node>(nid);
            std::vector<TCellID> adj_n = m_n2n[nid];
            if (adj_n.size() == 8) {
                // We've got 4 adjacent quad faces so
                // We apply an equipotential smoothing
                Node N[8];
                std::vector<Face> adj_faces = ni.get<Face>();
                std::vector<Node> f0_nodes = adj_faces[0].get<Node>();
                int position_ni=0;
                for(auto i=0; i<4;i++){
                    if(f0_nodes[i].id()==nid)
                        position_ni=i;
                }
                N[7] = f0_nodes[(position_ni+1)%4];
                N[0] = f0_nodes[(position_ni+2)%4];
                N[1] = f0_nodes[(position_ni+3)%4];

                //We look for the second faces sharing edge i-1 with f0
                bool found_f1=false;
                Face f1;
                for(auto i_f = 1; i_f<adj_faces.size() && !found_f1; i_f++){
                    Face fi = adj_faces[i_f];
                    std::vector<TCellID> fi_nids = fi.getIDs<Node>();
                    for(auto n_id:fi_nids){
                        if(n_id==N[1].id()){
                            found_f1 = true;
                            f1 = fi;
                        }
                    }
                }

                std::vector<Node> f1_nodes=f1.get<Node>();;
                position_ni=0;
                for(auto i=0; i<4;i++) {
                    if (f1_nodes[i].id() == nid)
                        position_ni = i;
                }
                if(f1_nodes[(position_ni+1)%4].id()==N[1].id()){
                    N[2] = f1_nodes[(position_ni+2)%4];
                    N[3] = f1_nodes[(position_ni+3)%4];
                }
                else{
                    N[2] = f1_nodes[(position_ni+2)%4];
                    N[3] = f1_nodes[(position_ni+1)%4];
                }

                //We look for the third faces sharing edge i-3 with f1
                bool found_f2=false;
                Face f2;
                for(auto i_f = 1; i_f<adj_faces.size() && !found_f2; i_f++){
                    Face fi = adj_faces[i_f];
                    if(fi!=f1) {
                        std::vector<TCellID> fi_nids = fi.getIDs<Node>();
                        for (auto n_id: fi_nids) {
                            if (n_id == N[3].id()) {
                                found_f2 = true;
                                f2=fi;
                            }
                        }
                    }
                }
                std::vector<Node> f2_nodes =f2.get<Node>();
                position_ni=0;
                for(auto i=0; i<4;i++) {
                    if (f2_nodes[i].id() == nid)
                        position_ni = i;
                }
                if(f2_nodes[(position_ni+1)%4].id()==N[3].id()){
                    N[4] = f2_nodes[(position_ni+2)%4];
                    N[5] = f2_nodes[(position_ni+3)%4];
                }
                else{
                    N[4] = f0_nodes[(position_ni+2)%4];
                    N[5] = f0_nodes[(position_ni+1)%4];
                }

                //We look for the fourth faces sharing edge i-7 with f0
                bool found_f3=false;
                Face f3;
                for(auto i_f = 1; i_f<adj_faces.size() && !found_f2; i_f++){
                    Face fi = adj_faces[i_f];
                    if(fi!=f1 && fi!=f2) {
                        found_f3 = true;
                        f3=fi;
                    }
                }
                std::vector<Node> f3_nodes =f3.get<Node>();
                position_ni=0;
                for(auto i=0; i<4;i++) {
                    if (f3_nodes[i].id() != nid &&
                        f3_nodes[i].id() != N[5].id() &&
                        f3_nodes[i].id() != N[7].id() ) {
                        N[6] = f0_nodes[(position_ni + 1) % 4];
                    }
                }


                math::Point P[8];
                for (auto i = 0; i < 8; i++) {
                    P[i] =N[i].point();
                }
                double xp = 0.5 * (P[1].X() - P[5].X());
                double yp = 0.5 * (P[1].Y() - P[5].Y());
                double zp = 0.5 * (P[1].Z() - P[5].Z());

                double xq = 0.5 * (P[7].X() - P[3].X());
                double yq = 0.5 * (P[7].Y() - P[3].Y());
                double zq = 0.5 * (P[7].Z() - P[3].Z());

                double alpha = xp * xp + yp * yp + zp * zp;
                double beta = xp * xq + yp * yq + zp * zq;
                double gamma = xq * xq + yq * yq + zq * zq;

                double W[8] = {-beta / 2.0,//W0
                               alpha,    //W1
                               beta / 2.0, //W2
                               gamma,    //W3
                               -beta / 2.0,//W4
                               alpha,    //W5
                               beta / 2.0, //W6
                               gamma}; //W7
                math::Point pi = ni.point();
                math::Point new_pi = pi;
                for (auto i = 0; i < 8; i++) {
                    new_pi = new_pi + W[i] * (P[i] - pi);
                }
                s->project(new_pi);
                ni.setPoint(new_pi);
            }
            else
            { //general case simple laplacian
                math::Point pi(0, 0, 0);
                for (auto adj_id: adj_n) {
                    math::Point pj = m->get<Node>(adj_id).point();
                    pi = pi + pj;
                }
                pi = pi * (1.0 / adj_n.size());
                try {
                    s->project(pi);
                }
                catch (GMDSException &e) { ; }
                ni.setPoint(pi);
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
