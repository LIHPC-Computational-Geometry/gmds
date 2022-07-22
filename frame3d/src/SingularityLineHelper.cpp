
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/math/Constants.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Triangle.h>
#include <gmds/frame3d/SingularityLineHelper.h>
#include <gmds/math/Line.h>
#include <gmds/math/Orientation.h>
/*---------------------------------------------------------------------------*/
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::frame3d;
/*---------------------------------------------------------------------------*/
SingularityLineHelper::SingularityLineHelper(Mesh* AMesh,
                                             ParamsGlobal& AGlobal,
                                             ParamsMark& AMarks)
        : m_mesh(AMesh)
        , m_params_gl(AGlobal)
        , m_bm(AMarks)
{
    m_chart_field = m_mesh->getVariable<math::Chart,GMDS_NODE>("SHChart");
    m_new_constraint = m_mesh->newVariable<math::Chart,GMDS_NODE>("chart_constraint");
    m_has_constraint = m_mesh->newVariable<int,GMDS_NODE>("has_constraint");
    for(auto n_id:m_mesh->nodes()){
        m_has_constraint->set(n_id,0);
    }
    math::Orientation::initialize();
}
/*---------------------------------------------------------------------------*/
SingularityLineHelper::~SingularityLineHelper()
{
    math::Orientation::finalize();

}
/*---------------------------------------------------------------------------*/
void
SingularityLineHelper::execute()
{
    //       initializeBooleanMarks();

    //=========================================================================
    // DETECTION OF SINGULAR ELEMENTS
    //=========================================================================
    detectSingularCells();

    Mesh lines(MeshModel(DIM3|F|N|F2N));
    for (auto f_id : m_sing_bnd_tri){
        std::cout<<"Detected boundary sing triangle: "<<f_id<<std::endl;
        Face f = m_mesh->get<Face>(f_id);
        m_sing_bnd_seps[f_id] = defineBoundarySlotsViaVectors(f, f.center());
        Node origin = lines.newNode(f.center());
        for(auto s:m_sing_bnd_seps[f_id]){
            Node to = lines.newNode(origin.point() + s);
            lines.newTriangle(origin,to,to);
        }
    }



    IGMeshIOService ioService(&lines);
    VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::F);
    writer.write("sweep_ff_separatrices.vtk");

    transportBndSingularities();


}
/*---------------------------------------------------------------------------*/
void
SingularityLineHelper::initializeBooleanMarks()
{
//        m_mark_line_sing = m_mesh->getNewMark<Region>();
//        m_mark_volume_pnt_sing = m_mesh->getNewMark<Region>();
//        m_mark_face_sing = m_mesh->getNewMark<Face>();
//        m_mark_edge_sing = m_mesh->getNewMark<Edge>();
//        m_mark_node_sing = m_mesh->getNewMark<Node>();
}
/*---------------------------------------------------------------------------*/
void
SingularityLineHelper::finalizeBooleanMarks()
{
//        m_mesh->unmarkAll<Region>(m_mark_volume_pnt_sing);
//        m_mesh->freeMark<Region>(m_mark_volume_pnt_sing);
//
//        m_mesh->unmarkAll<Region>(m_mark_line_sing);
//        m_mesh->freeMark<Region>(m_mark_line_sing);
//
//        m_mesh->unmarkAll<Face>(m_mark_face_sing);
//        m_mesh->freeMark<Face>(m_mark_face_sing);
//
//        m_mesh->unmarkAll<Edge>(m_mark_edge_sing);
//        m_mesh->freeMark<Edge>(m_mark_edge_sing);
//
//        m_mesh->unmarkAll<Node>(m_mark_node_sing);
//        m_mesh->freeMark<Node>(m_mark_node_sing);
}
/*----------------------------------------------------------------------------*/
void
SingularityLineHelper::detectSingularCells()
{
    for (auto r_id : m_mesh->regions()){
        Region current_region = m_mesh->get<Region>(r_id);
        std::vector<TCellID > node_ids = current_region.getIDs<Node>();
        bool onPnt = false;
        for (auto nid: node_ids){
            if (m_mesh->isMarked<Node>(nid, m_bm.mark_node_on_pnt))
                onPnt = true;
        }
        if (!onPnt) {
            std::vector<TCellID> node_ids = current_region.getIDs<Node>();
            math::Chart q[4];
            for(int i_n=0; i_n<4; i_n++){
                q[i_n]= (*m_chart_field)[node_ids[i_n]];
                }

            int sing_type = math::Chart::testSingularity(q[0],q[1],q[2],q[3]);
            if (sing_type != 0){
                m_sing_tet.push_back(current_region.id());
                std::vector<Face> fs = current_region.get<Face>();
                for (auto f : fs) {

                    std::vector<TCellID> node_ids = f.getIDs<Node>();
                    math::Chart qf[3];
                    for(int i_n=0; i_n<3; i_n++){
                        qf[i_n]=(*m_chart_field)[node_ids[i_n]];
                    }
                    int sing_type_f = math::Chart::testSingularity(qf[0],qf[1],qf[2]);
                    if (sing_type_f != 0) {
                        if (m_mesh->isMarked(f, m_bm.mark_face_on_surf)) {
                            m_sing_bnd_tri.push_back(f.id());

                            std::cout << "bnd type: " << computeSingularityIndex(f) << std::endl;


                        }
                    }
                }
            }
            m_sing_tet_type[r_id]=sing_type;

        }
    }
    std::cout << "Singular tets    : " << m_sing_tet.size() << ")\n";
    std::cout << "Singular bnd tris: " << m_sing_bnd_tri.size() << ")\n"<<std::endl;
}
/*----------------------------------------------------------------------------*/
int SingularityLineHelper::computeSingularityIndex(Face& AF)
{

    math::Vector3d normal = AF.normal();

    std::vector<Node> nodes = AF.get<Node>();

    std::vector<math::Chart> charts;
    charts.push_back(m_chart_field->value(nodes[0].id()));
    charts.push_back(m_chart_field->value(nodes[1].id()));
    charts.push_back(m_chart_field->value(nodes[2].id()));

    //We rotate the charts to be aligned along the face normal
    charts[0].align(normal);
    charts[1].align(normal);
    charts[2].align(normal);

    //Now we compute a reference vector angle for each node in
    //the face plane.

    std::vector<double> angle;
    angle.resize(3);
    std::vector<math::Vector3d> ref_vectors;
    ref_vectors.resize(3);
    //We take the first edge of the triangle as the OX axis
    math::Vector3d OX =nodes[1].point()- nodes[0].point();
    OX.normalize();

    for(auto i = 0; i < 3; i++) {
        //we select a chart vector which is not aligned with the
        //face normal
        double val_align = charts[i][0].dot(normal);
        int index_align = 0;
        for(auto j=1;j<3;j++){
            double val_j = charts[i][j].dot(normal);
            if(val_j>val_align){
                val_align = val_j;
                index_align=j;
            }
        }
        math::Vector3d v = charts[i][(index_align+1)%3];

        double theta = v.orientedAngle(OX,normal);

        angle[i] = math::modulo2PI(4*theta);
        math::AxisAngleRotation rotation(normal,angle[i]);
        ref_vectors[i]=rotation*OX;
    }

    math::Vector3d vi = ref_vectors[0];
    math::Vector3d vj = ref_vectors[1];
    math::Vector3d vk = ref_vectors[2];


    double wij = vi.orientedAngle(vj);
    double wjk = vj.orientedAngle(vk);
    double wki = vk.orientedAngle(vi);

    return std::round((wij+wjk+wki)/math::Constants::PI2);
}


/*----------------------------------------------------------------------------*/
void SingularityLineHelper::
updateSmoothingData(std::set<TCellID>& ANodesToUpdate,
                    Variable<math::Chart>* ANodeChartInfo,
                    Variable<int>* ANodeVarConstraint,
                    Variable<int>* ARegionsToUpdate)
{
    ANodesToUpdate.clear();
    int nb_reg=0;
    for (auto t_id:m_mesh->regions()) {
        Region t = m_mesh->get<Region>(t_id);
        (*ARegionsToUpdate)[t_id] = 0;
        std::vector<TCellID> t_ids = t.getIDs<Node>();
        if(isSingularTet(t_id)){
            if ((*ANodeVarConstraint)[t_ids[0]] != 1 ||
                (*ANodeVarConstraint)[t_ids[1]] != 1 ||
                (*ANodeVarConstraint)[t_ids[2]] != 1 ||
                (*ANodeVarConstraint)[t_ids[3]] != 1) {
                (*ARegionsToUpdate)[t_id] = 1;
                nb_reg++;
                ANodesToUpdate.insert(t_ids.begin(), t_ids.end());
            }
        }
    }
    std::cout<<"#### NB REG TO UPDATE: "<<nb_reg<<std::endl;
}

/*----------------------------------------------------------------------------*/
void SingularityLineHelper::
initFrameSmoothing(std::set<TCellID>& ANodesToUpdate,
                   Variable<int>* ANodeVarConstraint,
                   Variable<int>* ARegionsToUpdate)
{
    std::set<TCellID> constraint_nodes;
    for(auto n_id:m_mesh->nodes()){
        if((*ANodeVarConstraint)[n_id]==1){
            constraint_nodes.insert(n_id);
        }
    }



    for(auto r_id : m_mesh->regions()){
        if((*ARegionsToUpdate)[r_id]==0)
            continue;

        TCellID bnd_face_id[2];

        Region current_tet = m_mesh->get<Region>(r_id);
        std::vector<Node> current_nodes = current_tet.get<Node>();

        math::Chart ci = (*m_chart_field)[current_nodes[0].id()];


        double d0 = closestBndFace(r_id,  ci.Z(), bnd_face_id[0]);
        double d1 = closestBndFace(r_id, -ci.Z(), bnd_face_id[1]);

        TCellID  closest_bnd_face_id=NullID;
        if(d0<d1){
            closest_bnd_face_id=bnd_face_id[0];
        }
        else{
            closest_bnd_face_id=bnd_face_id[1];
        }

        Face bnd_face = m_mesh->get<Face>(closest_bnd_face_id);
        std::vector<TCellID> bnd_nodes = bnd_face.getIDs<Node>();
        std::vector<math::Quaternion>  q;
        std::vector<double> w;
        for(auto bnd_nid:bnd_nodes){
            q.push_back(math::Quaternion((*m_chart_field)[bnd_nid]));
            w.push_back(1);
        }
//        math::Chart bnd_chart(math::Quaternion::mean(q,w));
        math::Chart bnd_chart = (*m_chart_field)[bnd_nodes[1]];
        for(auto n:current_nodes){
            math::Chart cn = (*m_chart_field)[n.id()];
            math::Chart bnd_copy = bnd_chart;
          //  bnd_copy.align(cn.Z());
            (*m_chart_field)[n.id()]= bnd_copy;
        }

    }
}
/*----------------------------------------------------------------------------*/
bool SingularityLineHelper::
isSingularTet(const TCellID AId)
{
    Region t = m_mesh->get<Region>(AId);
    std::vector<TCellID> t_ids = t.getIDs<Node>();
    math::Chart c0 = (*m_chart_field)[t_ids[0]];
    math::Chart c1 = (*m_chart_field)[t_ids[1]];
    math::Chart c2 = (*m_chart_field)[t_ids[2]];
    math::Chart c3 = (*m_chart_field)[t_ids[3]];
    return (math::Chart::testSingularity(c0, c1, c2, c3) != 0);

}
/*----------------------------------------------------------------------------*/
void SingularityLineHelper::transportBndSingularities()
{
    transportBndSingularity(m_sing_bnd_tri[1]);

    return;
    for (auto f_id : m_sing_bnd_tri){
  /*      if(f_id==65527 || f_id==78756){
            continue;
        }
        if(f_id==57908){
            moveBoundarySingularity(f_id,61727);
            f_id = 61727;
        }
    */    transportBndSingularity(f_id);
    }

}
/*----------------------------------------------------------------------------*/
void SingularityLineHelper::
moveBoundarySingularity(TCellID AFrom, TCellID ATo)
{

    std::cout<<"Move singularity from face "<<AFrom<<std::endl;
    Face from_face = m_mesh->get<Face>(AFrom);
    std::cout<<"to face "<<AFrom<<std::endl;
    Face to_face = m_mesh->get<Face>(ATo);
    std::vector<Node> to_nodes = to_face.get<Node>();

    // we build an orthonormal basis local to to_face.
    math::Vector3d v1=to_nodes[1].point()- to_nodes[0].point();
    Region to_r = to_face.get<Region>()[0];
    math::Vector3d v3 = getInputNormal(to_face,to_r);
    math::Vector3d v2 = v1.cross(v3);
    math::Chart from_chart(v1,v2,v3);

    // convert separatrices into angle values in (v1,v2,v3) base
    std::vector<math::Vector3d> from_seps = m_sing_bnd_seps[from_face.id()];
    std::vector<double> sep_angles;
    sep_angles.reserve(from_seps.size());
    for(auto sep:from_seps){
        sep_angles.push_back(from_chart.X().orientedAngle(sep,from_chart.Z()));
    }


    m_sing_bnd_seps[to_face.id()]=from_seps;

    //By this way, we define new values at the corners of to_face
    for(auto ni:to_nodes){
        //We compute for node ni, its angle into the local chart
        math::Point pi = ni.point();
        math::Point ref_point = to_face.center();


        math::Vector3d n_vec=pi-ref_point;
        double n_angle = from_chart.X().orientedAngle(n_vec,from_chart.Z());

        bool found_interval=false;
        int index_min=0, index_max = 0;
        for(unsigned int i_angle=0; i_angle<sep_angles.size() && !found_interval; i_angle++){
            double val_i = sep_angles[i_angle];
            double val_j = sep_angles[(i_angle+1)%sep_angles.size()];
            if(val_i<=n_angle && n_angle<val_j){
                index_min=i_angle;
                index_max=(i_angle+1)%sep_angles.size();
                found_interval=true;
            }
        }
        if(!found_interval){
            //   std::cout<<"\t not found!!!"<<std::endl;
            //means smaller than the smallest one
            //or greater than the greatest one.
            index_min=sep_angles.size()-1;
            index_max=0;
        }

        double min_angle = sep_angles[index_min];

        math::Vector3d vmin = math::AxisAngleRotation(from_chart.Z(),min_angle)*from_chart.X();
        double max_angle = sep_angles[index_max];
        math::Vector3d vmax = math::AxisAngleRotation(from_chart.Z(),max_angle)*from_chart.X();
        vmin.normalize();
        vmax.normalize();

        double ref_angle = vmin.orientedAngle(vmax,from_chart.Z());
        double val_angle = vmin.orientedAngle(n_vec,from_chart.Z());

        double ratio = std::abs(val_angle/ref_angle);

        vmin = vmin.cross(from_chart.Z());
        if(vmax.dot(vmin)<0){
            vmin= vmin.opp();
        }
        vmin.normalize();
        vmax.normalize();
        math::Vector3d v = (1-ratio)*vmin + ratio*vmax;
        v.normalize();
        math::Vector3d v2  = v.cross(from_chart.Z());
        v2.normalize();

        //We update the chart attached to the AToNodes
        (*m_chart_field)[ni.id()]=math::Chart(v,v2,from_chart.Z());
    }


}
/*----------------------------------------------------------------------------*/
void SingularityLineHelper::transportBndSingularity(gmds::TCellID &AFID)
{
    Face f = m_mesh->get<Face>(AFID);
    std::cout<<"Start from face "<<AFID<<std::endl;
    std::vector<Node> fns = f.get<Node>();
   // std::vector<math::Chart> node_charts;

    //only one region along the boundary
    Region r = f.get<Region>()[0];

    math::Vector3d v = getInputNormal(f,r);

    PointVolumetricData pvd(f.center(),v,r);

    math::Point p;
    std::vector<TCellID> tets;
    std::map<TCellID,math::Chart> singularities;
    followFlow(f,pvd,p,tets, singularities);

    Variable<math::Vector3d>* vx = NULL;
    Variable<math::Vector3d>* vy = NULL;
    Variable<math::Vector3d>* vz = NULL;

    try{
        vx = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("new_sing_X");
        vy = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("new_sing_Y");
        vz = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("new_sing_Z");
    }
    catch (GMDSException&e){
        vx = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("new_sing_X");
        vy = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("new_sing_Y");
        vz = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("new_sing_Z");

    }


    for(auto sing:singularities){
        TCellID n_id = sing.first;
        Node n= m_mesh->get<Node>(n_id);
        std::vector<Region> n_tets = n.get<Region>();
        vx->set(sing.first,sing.second.X());
        vy->set(sing.first,sing.second.Y());
        vz->set(sing.first,sing.second.Z());
        m_has_constraint->set(sing.first,1);
        m_new_constraint->set(sing.first,sing.second);
    }

    Variable<int>* v_line  = m_mesh->newVariable<int, GMDS_REGION>("line_"+std::to_string(AFID));
    Variable<int>* v_sing  = m_mesh->newVariable<int, GMDS_REGION>("line_sing"+std::to_string(AFID));

    for(auto i:tets){
        v_line->set(i,1);

    }
    int sing_tets=0;
    for(auto i:tets) {

        std::vector<TCellID> node_ids = m_mesh->get<Region>(i).getIDs<Node>();
        math::Chart q[4];
        for (int i_n = 0; i_n < 4; i_n++) {
            if((*m_has_constraint)[node_ids[i_n]]==0){
                std::cout<<"\t error! No constraint on node "<<node_ids[i_n]<<" in tet "<<i<<std::endl;
            }
            q[i_n] = (*m_new_constraint)[node_ids[i_n]];
        }

        int sing_type = math::Chart::testSingularity(q[0], q[1], q[2], q[3]);
        if(sing_type!=0){
            sing_tets++;
        }
        v_sing->set(i,sing_type);
    }
    std::cout<<"NB (SING/TRAVERS. TETS) = "<<sing_tets<<" / "<<tets.size()<<std::endl;
    std::cout<<"reach point "<<p<<std::endl;
}
/*----------------------------------------------------------------------------*/
double SingularityLineHelper::
closestBndFace(const TCellID ATetID, math::Vector3d AV, TCellID& ABndFaceID)
{
    double dist=0;
    Region current_cell = m_mesh->get<Region>(ATetID);
    math::Vector3d from_direction = AV;
    math::Point from_point = current_cell.center();

    TCellID from_face_id=NullID;
    math::Point prev_point = from_point;
    bool reach_bnd=false;
    while (reach_bnd==false) {
        //to avoid starting in the wrong being on an edge
        from_point  = 0.95*from_point+ 0.05*current_cell.center();

        //================================================
        // We build the triangular shapes corresponding to
        // each tet face
        std::vector<Face> cf = current_cell.get<Face>();

        std::vector<math::Triangle> t;
        t.resize(4);
        for(unsigned int i=0; i<4; i++){
            std::vector<Node> cn = cf[i].get<Node>();
            t[i]= math::Triangle(cn[0].point(),
                                 cn[1].point(),
                                 cn[2].point());
        }
        //================================================
        // PHASE 1
        //================================================
        math::Point out_pnt;
        math::Vector3d out_dir;
        int out_index=-1;
        computeFuzzyHeuns(from_point, from_direction, cf, t,
                          out_pnt, out_dir, out_index);
        //================================================
        // PHASE 2
        //================================================
        from_direction = 0.5*from_direction + 0.5*out_dir;
        computeFuzzyHeuns(from_point, from_direction, cf, t,
                          out_pnt, out_dir, out_index);

        //================================================
        // FINALIZATION
        //================================================
        // Now we have our point to get out of current_cell
        math::Vector3d v= out_pnt-from_point;

        Face out_face = cf[out_index];

        //Fuzzy approach to avoid topological cases
        math::Point new_from =0.95*out_pnt+ 0.05*out_face.center();

        //we update the direction we come from for the next step
        from_direction = out_dir;
        dist +=(new_from-from_point).norm();

        from_point  = new_from;


        //We compute the next current cell using the face we execute
        //from_point

        std::vector<Region> adj_out_face = out_face.get<Region>();
        if(adj_out_face.size()==1){
            //we are on the boundary even if we have not reach the
            //expected distance
            reach_bnd=true;
            ABndFaceID = out_face.id();
        }
        else{
            if(adj_out_face[0].id()==current_cell.id())
                current_cell=adj_out_face[1];
            else
                current_cell=adj_out_face[0];
        }
        //As we arrive by a face, we only have one output point
        from_face_id=out_face.id();

    };


    return dist;
}
/*----------------------------------------------------------------------------*/
bool SingularityLineHelper::followFlow(const Face AFromFace,
                                       const PointVolumetricData& AData,
                                       math::Point&               APnt,
                                       std::vector<TCellID>& ACrossedTet,
                                       std::map<TCellID,math::Chart>& ASingVal)
{
    Face from_face = AFromFace;


    /*Each node of an intersected tet keep in mind:
     * - the propagation direction, i.e the orthogonal plane normal it belongs to,
     * - a chart defining a right-hand sided basis for transporting singularity,
     * - an intersection point between the plane it belongs to and the propag direction.
    */
    Variable<math::Vector3d>* propag_direction=NULL;
    Variable<math::Chart>*    propag_chart =NULL;
    Variable<math::Point>*    propag_intersection=NULL;

    int mark_node_assigned = m_mesh->newMark<Node>();
    try{
        propag_direction    = m_mesh->newVariable<math::Vector3d, GMDS_NODE>("propag_direction");
        propag_chart        = m_mesh->newVariable<math::Chart   , GMDS_NODE>("propag_chart");
        propag_intersection = m_mesh->newVariable<math::Point   , GMDS_NODE>("propag_intersection");

    }
    catch(GMDSException &e){
        propag_direction    = m_mesh->getVariable<math::Vector3d, GMDS_NODE>("propag_direction");
        propag_chart        = m_mesh->getVariable<math::Chart   , GMDS_NODE>("propag_chart");
        propag_intersection = m_mesh->getVariable<math::Point   , GMDS_NODE>("propag_intersection");

    }
    // We build the reference chart to propagate singularity line directions
    std::vector<Node> from_nodes = from_face.get<Node>();
    math::Vector3d v1=from_nodes[1].point()- from_nodes[0].point();
    math::Vector3d v3(AData.dir);
    math::Vector3d v2 = v1.cross(v3);
    math::Chart from_chart(v1,v2,v3);

    Mesh lines(MeshModel(DIM3|F|N|F2N));
    Mesh path(MeshModel(DIM3|F|N|F2N));


    //convert separatrices into angle values in (v1,v2) base
    std::vector<math::Vector3d> from_seps = m_sing_bnd_seps[from_face.id()];
    std::vector<double> sep_angles;
    sep_angles.reserve(from_seps.size());
    for(auto sep:from_seps){

        sep_angles.push_back(from_chart.X().orientedAngle(sep,from_chart.Z()));
    }

    std::sort(sep_angles.begin(),sep_angles.end());

    for(auto n:from_nodes){
        propag_chart->set(n.id(),from_chart);
        propag_intersection->set(n.id(),AData.pnt);
    }

    ACrossedTet.clear();
    double current_dist=0;
    math::Vector3d from_direction(AData.dir);
    math::Point from_point = AData.pnt;
    gmds::Region current_cell = AData.tet;

    math::Point prev_point = from_point;
    bool reach_bnd=false;
    while (reach_bnd==false) {
        ACrossedTet.push_back(current_cell.id());
        //to avoid starting in the wrong being on an edge
        from_point  = 0.95*from_point+ 0.05*current_cell.center();

        //================================================
        // We build the triangular shapes corresponding to
        // each tet face
        std::vector<Face> cf = current_cell.get<Face>();

        std::vector<math::Triangle> t;
        t.resize(4);
        for(unsigned int i=0; i<4; i++){
            std::vector<Node> cn = cf[i].get<Node>();
            t[i]= math::Triangle(cn[0].point(),
                                 cn[1].point(),
                                 cn[2].point());

            if(cf[i].id()==from_face.id()){
                for(auto n:cn){
                    (*propag_direction)[n.id()]=from_direction;

                    math::Plane plane(n.point(), from_direction);
                    (*propag_intersection)[n.id()]=plane.project(from_point);
                    Node origin = lines.newNode((*propag_intersection)[n.id()]);

                    math::Chart to_chart = from_chart;
                    to_chart.align(from_direction);
                    (*propag_chart)[n.id()]=to_chart;



                    for(auto angle:sep_angles){
                        math::Vector3d sep = math::AxisAngleRotation(to_chart.Z(),angle)*to_chart.X();
                        Node to = lines.newNode(origin.point() + sep);
                        lines.newTriangle(origin,to,to);
                    }
                    m_mesh->mark(n,mark_node_assigned);
                }


            }
        }

        //================================================
        // PHASE 1
        //================================================
        math::Point out_pnt;
        math::Vector3d out_dir;
        int out_index=-1;
        computeFuzzyHeuns(from_point, from_direction, cf, t,
                          out_pnt, out_dir, out_index);
        //================================================
        // PHASE 2
        //================================================
        from_direction = 0.5*from_direction + 0.5*out_dir;
        computeFuzzyHeuns(from_point, from_direction, cf, t,
                          out_pnt, out_dir, out_index);

        //================================================
        // FINALIZATION
        //================================================
        // Now we have our point to get out of current_cell
        math::Vector3d v=out_pnt-from_point;

        Face out_face = cf[out_index];
        //       std::cout<<"Out face: "<<out_face.getID()<<std::endl;
        //Fuzzy approach to avoid topological cases
        math::Point new_from =0.95*out_pnt+ 0.05*out_face.center();

        //we update the direction we come from for the next step
        from_direction = out_dir;
        current_dist +=(new_from-from_point).norm();

        Node n_from = path.newNode(from_point);
        Node n_to = path.newNode(new_from);
        path.newTriangle(n_from,n_to,n_to);
        from_point  = new_from;


        //We compute the next current cell using the face we execute
        //from_point

                std::vector<Region> adj_out_face = out_face.get<Region>();
        if(adj_out_face.size()==1){
            //we are on the boundary even if we have not reach the
            //expected distance
            APnt = from_point;
            reach_bnd=true;
        }
        else{
            if(adj_out_face[0].id()==current_cell.id())
                current_cell=adj_out_face[1];
            else
                current_cell=adj_out_face[0];
        }
        //As we arrive by a face, we only have one output point
        from_face=out_face;

        //and we recompute the chart
        from_chart.align(from_direction);
    };




    IGMeshIOService ioService(&lines);
    VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::F);
    writer.write("sweep_ff_generated_separatrices_"+std::to_string(AFromFace.id())+".vtk");

    IGMeshIOService ioServicep(&path);
    VTKWriter writerp(&ioServicep);
    writerp.setCellOptions(gmds::N|gmds::F);
    writerp.write("sweep_ff_generated_path_"+std::to_string(AFromFace.id())+".vtk");


    /*
      * Now we update frames along the path of tets we crossed.
      * By construction, they are ordered in the traversal direction
      */
    int mark_done = m_mesh->newMark<Node>();
    std::vector<TCellID> vec_done;


    // We keep in mind the 3 nodes of the face we come from to get into the current
    // tet
    std::vector<Node> incoming_nodes = AFromFace.get<Node>();


    for(unsigned int i_tet=0; i_tet<ACrossedTet.size();i_tet++) {
        Region tet = m_mesh->get<Region>(ACrossedTet[i_tet]);
        std::vector<Node> tet_nodes = tet.get<Node>();
        for (auto n:tet_nodes) {
            if (!m_mesh->isMarked(n, mark_node_assigned)) {
                std::cout << "\t ERROR NODE " << n.id() << " in tet " << tet.id() << std::endl;
                continue;
            }
        }
        //We look for the incoming nodes that were already assigned and we only modify the fourth one
        Node current_node;
        for(auto n_tet:tet_nodes){
            bool found_node = false;
            for(auto n_face:incoming_nodes){
                if(n_tet.id()==n_face.id()){
                    found_node=true;
                }
            }
            if(!found_node){
                current_node=n_tet;
            }
        }


        //We only modify the current node
        math::Vector3d plane_direction = (*propag_direction)[current_node.id()];
        math::Chart plane_chart     = (*propag_chart)[current_node.id()];
        math::Point plane_point     = (*propag_intersection)[current_node.id()];


        math::Vector3d n_vec=current_node.point()-plane_point;
        double n_angle = plane_chart.X().orientedAngle(n_vec,plane_chart.Z());

        bool found_interval=false;
        int index_min=0, index_max = 0;
        for(unsigned int i_angle=0; i_angle<sep_angles.size() && !found_interval; i_angle++){
            double val_i = sep_angles[i_angle];
            double val_j = sep_angles[(i_angle+1)%sep_angles.size()];
            if(val_i<=n_angle && n_angle<val_j){
                index_min=i_angle;
                index_max=(i_angle+1)%sep_angles.size();
                found_interval=true;
            }
        }
        if(!found_interval){
         //   std::cout<<"\t not found!!!"<<std::endl;
            //means smaller than the smallest one
            //or greater than the greatest one.
            index_min=sep_angles.size()-1;
            index_max=0;
        }

        double min_angle = sep_angles[index_min];

        math::Vector3d vmin = math::AxisAngleRotation(plane_chart.Z(),min_angle)*plane_chart.X();
        double max_angle = sep_angles[index_max];
        math::Vector3d vmax = math::AxisAngleRotation(plane_chart.Z(),max_angle)*plane_chart.X();
        vmin.normalize();
        vmax.normalize();

        double ref_angle = vmin.orientedAngle(vmax,from_chart.Z());
        double val_angle = vmin.orientedAngle(n_vec,from_chart.Z());

        double ratio = std::abs(val_angle/ref_angle);

        vmin = vmin.cross(plane_chart.Z());
        if(vmax.dot(vmin)<0){
            vmin= vmin.opp();
        }
        vmin.normalize();
        vmax.normalize();
        math::Vector3d v = (1-ratio)*vmin + ratio*vmax;
        v.normalize();
        math::Vector3d v2  = v.cross(plane_chart.Z());
        v2.normalize();


        ASingVal[current_node.id()]=math::Chart(v,v2,plane_chart.Z());
    //    std::cout<<current_node.id()<<") "<<v<<", "<<v2<<", "<<plane_chart.Z()<<std::endl;
        //We build the list of nodes shared by the current face and the next one
        if(i_tet!=ACrossedTet.size()-1){
            Region next_tet = m_mesh->get<Region>(ACrossedTet[i_tet+1]);
            std::vector<TCellID > next_node_ids= next_tet.getIDs<Node>();
            incoming_nodes.clear();
            for(auto n_tet:tet_nodes){
                bool found_node = false;
                for(auto n_id:next_node_ids){
                    if(n_tet.id()==n_id){
                        found_node = true;
                    }
                }
                if(found_node){
                    incoming_nodes.push_back(n_tet);
                }
            }
        }
    }

    std::vector<TCellID>  init_node_ids = AFromFace.getIDs<Node>();
    for(auto i:init_node_ids){
        ASingVal[i]=(*m_chart_field)[i];
    }
    for(auto i:vec_done){
        m_mesh->unmark<Node>(i,mark_done);
    }
    m_mesh->freeMark<Node>(mark_done);
    for(auto tet_id:ACrossedTet){
        Region tet = m_mesh->get<Region>(tet_id);
        std::vector<TCellID > tet_node_ids = tet.getIDs<Node>();
        for(auto n_id:tet_node_ids){
            m_mesh->unmark<Node>(n_id, mark_node_assigned);
        }
    }
    m_mesh->freeMark<Node>(mark_node_assigned);

    return true;
}

/*----------------------------------------------------------------------------*/
void SingularityLineHelper::
computeFuzzyHeuns(const math::Point&                 AFromPnt,
                  const math::Vector3d&              AFromDir,
                  const std::vector<Face>&           AFaces,
                  const std::vector<math::Triangle>& ATri,
                  math::Point&                       AToPnt,
                  math::Vector3d&                    AToDir,
                  int&                               AToFaceId)
{

    // check whether the line intersects the triangle
    math::Line ray(AFromPnt,math::Vector3d({AFromDir.X(), AFromDir.Y(), AFromDir.Z()}));

    //===================================================================
    double param[4] = {-1, -1,-1,-1};
    math::Point p[4];

    for(auto i=0; i<4; i++){
        math::Plane pli = ATri[i].getPlaneIncluding();
        if(ray.intersect3D(pli, p[i], param[i])){
            //intersect the plane, but the triangle??
            //before computing the bar coordinate, we eliminate intersection
            //at the infinity almost due to almost parallel ray and plane
            std::vector<Node> n = AFaces[i].get<Node>();
            math::Point pn[3] ={
                    n[0].point(),
                    n[1].point(),
                    n[2].point()};
            //We compute the barycentric coords.
            bool on_edge[3]={false,false,false};
            if(!isIn(p[i],AFaces[i], on_edge[0],on_edge[1],on_edge[2])){
                param[i]=-1;

            };

        }
        else{
            param[i]=-1;
        }
//          std::cout<<"\t intersection with face "<<AFaces[i].getID()<<" -> "<<param[i]<<std::endl;
    }
    double best_param =param[0];
    auto out_index = 0;

    for(auto i=1; i<4; i++){
        if(param[i]>best_param){
            out_index=i;
            best_param=param[i];
        }
    }
    if(best_param<=1e-8)
        throw GMDSException("Tools::computeFuzzyHeuns: No out face (1)");
    AToPnt=p[out_index];
//    std::cout<<"\t OUT PNT: "<<AToPnt<<std::endl;
    //===================================================================
    //compute the frame in out_pnt
    math::Chart ci = computeChartIn(AToPnt, AFaces[out_index]);

    //===================================================================
    //among the 6 vectors of ci, we take the one which is the
    // best aligned with dir and we start the process a second
    // time
    math::Vector3d ci_vectors[6] = {
            ci.X(), -ci.X(),  ci.Y(),
            -ci.Y(),  ci.Z(), -ci.Z()
    };
    math::Vector3d heuns_corr = ci_vectors[0];
    double best_align_dot = AFromDir.dot(ci_vectors[0]);

    for(int i=0; i<6;i++){
        if(best_align_dot<AFromDir.dot(ci_vectors[i])){
            heuns_corr = ci_vectors[i];
            best_align_dot = AFromDir.dot(ci_vectors[i]);
        }
    }

    AToDir= heuns_corr;
    AToFaceId = out_index;
}
math::Chart SingularityLineHelper::computeChartIn(const math::Point& APnt,
                                                  const Face&        AFace)
{
    std::vector<Node> n =AFace.get<Node>();
    //===================================================================
    //STEP 1 - We compute the location of APnt into AFace
    //===================================================================
    double coeff[3]={0, 0, 0};

    math::Point::computeBarycentric(n[0].point(), n[1].point(),
                                    n[2].point(), APnt,
                                    coeff[0], coeff[1], coeff[2]);

    //===================================================================
    //STEP 2 - We extract the quaternion representation of each frame  at
    //         the face corners
    //===================================================================
    std::vector<math::AxisAngleRotation> r;
    r.resize(3);
    for(int i=0;i<3;i++){
        math::Chart ci = m_chart_field->value(n[i].id());
        r[i]=math::AxisAngleRotation(ci);
    }

    std::vector<math::Quaternion> qs;
    qs.resize(3);
    for(int i=0;i<3;i++){
        qs[i]=math::Quaternion(r[i].toChart());
    }

    std::vector<TCoord> ws;
    ws.resize(3);
    ws[0]=coeff[0];
    ws[1]=coeff[1];
    ws[2]=coeff[2];

    //===================================================================
    //STEP 3 - We compute the mean quaternion and return the corresponding
    //         chart
    //===================================================================
    math::Quaternion q = math::Quaternion::mean(qs,
                                                ws);

    return math::Chart(q);

}

/*---------------------------------------------------------------------------*/
bool SingularityLineHelper::isIn(const math::Point& AP,
                                 const Face& ATri,
                                 bool& AOnEdge0,
                                 bool& AOnEdge1,
                                 bool& AOnEdge2)
{
    std::vector<Node> n = ATri.get<Node>();

    math::Point pnt[3] = {
            n[0].point(),
            n[1].point(),
            n[2].point()
    };
    math::Plane plane(pnt[0],pnt[1],pnt[2]);
    math::Point p = plane.project(AP);

    // We look if p is insie or outside of the triangle defined by ATri
    math::Vector3d normal = plane.getNormal();
    normal.normalize();
    math::Orientation::Sign ori[3] ={
            math::Orientation::orient3d(p, pnt[1] , pnt[2] , p+normal),
            math::Orientation::orient3d(p, pnt[2] , pnt[0] , p+normal),
            math::Orientation::orient3d(p, pnt[0] , pnt[1] , p+normal)
    };
    if ((ori[0] >= 0 && ori[1] >= 0 && ori[2] >= 0 ) ||
        (ori[0] <= 0 && ori[1] <= 0 && ori[2] <= 0 ) ) {
        if(ori[0]==math::Orientation::ZERO){
            AOnEdge0=true;
        }
        if(ori[1]==math::Orientation::ZERO){
            AOnEdge1=true;
        }
        if(ori[2]==math::Orientation::ZERO){
            AOnEdge2=true;
        }
        return true;
    }
    return false;
}
/*----------------------------------------------------------------------------*/
math::Vector3d
SingularityLineHelper::getOutputNormal(const Face& AF, const Region& AR)
{
    std::vector<Node> r_nodes = AR.get<Node>();
    std::vector<Node> f_nodes = AF.get<Node>();

    if (r_nodes.size() != 4)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
    if (f_nodes.size() != 3)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");

    // we execute through all the nodes of ARegion to find the one that do not belong
    // to AF
    for (auto n : r_nodes) {
        if (n != f_nodes[0] && n != f_nodes[1] && n != f_nodes[2]) {
            // n is the node opposite to the face AF
            math::Vector3d normal = AF.normal();
            math::Vector3d in_vector= n.point()-f_nodes[0].point();
            if (normal.dot(in_vector) > 0.0) {
                return math::Vector3d({-normal.X(), -normal.Y(), -normal.Z()});
            } else {
                return normal;
            }
        }
    }
    throw GMDSException("SingularityGraphBuilder::getOutputNormal unexpected behaviour");
}
/*----------------------------------------------------------------------------*/
math::Vector3d
SingularityLineHelper::getInputNormal(const Face& AF, const Region& AR)
{
    math::Vector3d outVec = getOutputNormal(AF, AR);
    return outVec.opp();
}
/*----------------------------------------------------------------------------*/
std::vector<math::Vector3d>
SingularityLineHelper::
defineBoundarySlotsViaAngles(const Face& ABndFace,
                             const math::Point& ASingLoc)
{
    std::vector<math::Vector3d> sep;

    std::vector<TCellID> node_ids = ABndFace.getIDs<Node>();
    std::vector<Node>    nodes    = ABndFace.get<Node>();


    math::Vector3d normal = ABndFace.normal();
    //=====================================================================
    // We select on chart vector living in the plane of f
    math::Chart c[3] = {(*m_chart_field)[node_ids[0]],
                        (*m_chart_field)[node_ids[1]],
                        (*m_chart_field)[node_ids[2]]};
    math::Vector3d v[3];
    math::Vector3d n3d(normal);
    for (auto i = 0; i < 3; i++) {
        for (auto j = 0; j < 3; j++) {
            if (abs(c[i][j].dot(n3d)) < 0.1) {
                v[i] = c[i][j];
            }
        }
    }

    //=====================================================================
    // We compute the singularity index (+1 or -1)
    //=====================================================================
    // We compute reference angle in the face plane using Palacio technique
    math::Vector3d ref=nodes[1].point()-nodes[0].point();

    ref.normalize();

    double angle[3] = {0, 0, 0};
    for (auto i = 0; i < 3; i++) {
        TCoord a = v[i].orientedAngle(ref, normal);
        if (a < 0) {
            a = math::Constants::PI2 + a;
        }
        angle[i] = math::modulo2PI(4 * a);
    }

    //=====================================================================
    // Then the reference vector
    math::Vector3d ref0({cos(angle[0]), sin(angle[0]), 0.0});
    math::Vector3d ref1({cos(angle[1]), sin(angle[1]), 0.0});
    math::Vector3d ref2({cos(angle[2]), sin(angle[2]), 0.0});

    //=====================================================================
    // And so the face index k = 1 or -1 depending on the sing type
    double w01 = ref0.orientedAngle(ref1);
    double w12 = ref1.orientedAngle(ref2);
    double w20 = ref2.orientedAngle(ref0);
    double index_d = (w01 + w12 + w20) / math::Constants::PI2;
    int index = round(index_d);

    //=====================================================================
    // As we have the index, we can deduce singularity direction as done
    // in H. Fogg PHD mansucrit (see Bunin's papers too)

    // we need to find a first direction
    bool found_first_separatrix = false;
    math::Vector3d first_sepa;
    for (auto i = 0; i < 3 && !found_first_separatrix; i++) {
        const auto j = (i + 1) % 3;
        math::Point pi = nodes[i].point();
        math::Point pj = nodes[j].point();
        //        std::cout<<"From "<<nodes[i].getID()<<" to "<<nodes[j].getID()<<std::endl;
        // we work along the edge [i,j]
        math::Vector3d vij=pj-pi;

        // all ref angles are recomputed for this edge
        for (auto k = 0; k < 3; k++) {
            TCoord a = vij.orientedAngle(v[k], normal);
            if (a < 0) {
                a = math::Constants::PI2 + a;
            }
            angle[k] = math::modulo2PI(4 * a);
        }

        // We get the four angles of the cross in point i and j
        // cross in the plane of f

        double cross_angle_i[4];
        TCoord a = angle[i] / 4.0;
        for (auto k = 0; k < 4; k++) {
            cross_angle_i[k] = math::modulo2PI(a + k * math::Constants::PIDIV2);
        }

        double cross_angle_j[4];
        a = angle[j] / 4.0;
        for (auto k = 0; k < 4; k++) {
            cross_angle_j[k] = math::modulo2PI(a + k * math::Constants::PIDIV2);
        }
        // For each vector, we check if we find a point where to
        // execute out along edge [i,j]
        math::Vector3d vi=ASingLoc-pi;
        math::Vector3d vj=ASingLoc-pj;
        double alpha_i = vij.orientedAngle(vi, normal);
        double alpha_j = vij.orientedAngle(vj, normal);

        for (auto k = 0; k < 4 && !found_first_separatrix; k++) {
            double ai = cross_angle_i[k];
            double aj = cross_angle_j[0];

            double dij = abs(aj - ai);
            if (dij > math::Constants::PI)
                dij = math::Constants::PI2 - dij;
            int match_j = 0;
            for (auto l = 1; l < 4; l++) {
                double dij_l = abs(ai - cross_angle_j[l]);
                if (dij_l > math::Constants::PI)
                    dij_l = math::Constants::PI2 - dij_l;

                if (dij_l < dij) {
                    dij = dij_l;
                    aj = cross_angle_j[l];
                    match_j = l;
                }
            }
            double delta = aj - ai;
            double t = (alpha_i - ai) / (delta - (alpha_j - alpha_i));
            if (t >= 0 && t <= 1) {
                found_first_separatrix = true;
                math::Point pnt_sepa = (1 - t) * pi + t * pj;
                first_sepa = pnt_sepa-ASingLoc;
                first_sepa.normalize();
            }
        }

    }  // for(auto i=0; i<3 && !found_first_separatrix; i++)

    if (!found_first_separatrix) {
        throw GMDSException("Impossible to compute a separatrix");
    }

    double separatrix_angle = math::Constants::PI2 / (4 + index);

    int nb_sing = 4 + index;
    sep.resize(nb_sing);
    sep[0] = first_sepa;
    math::AxisAngleRotation rot_sepa(normal, separatrix_angle);
    for (auto s = 1; s < nb_sing; s++) {
        sep[s] = rot_sepa * sep[s - 1];
    }

    return sep;
}
/*----------------------------------------------------------------------------*/
std::vector<math::Vector3d>
SingularityLineHelper::
defineBoundarySlotsViaVectors(const Face& ABndFace,
                              const math::Point& ASingLoc)
{
    math::Vector3d normal = ABndFace.normal();
    std::vector<math::Vector3d> sep;

    double nx = normal.X();
    double ny = normal.Y();
    double nz = normal.Z();
    double cx = ASingLoc.X();
    double cy = ASingLoc.Y();
    double cz = ASingLoc.Z();

    //=====================================================================
    // We select on chart vector living closed to the plane of f
    math::Vector3d v[3];
    math::Vector3d n3d(normal);

    std::vector<Node> nodes = ABndFace.get<Node>();
    for (auto i = 0; i < nodes.size(); i++) {
        math::Chart c =(*m_chart_field)[nodes[i].id()];
        for (auto j = 0; j < 3; j++) {
            if (abs(c[j].dot(n3d)) < 0.1) {
                v[i] = c[j];
            }
        }
    }
    //=====================================================================
    // We project charts into the plane to avoid numerical issues
    math::Plane pl(ASingLoc, normal);
    for (auto i = 0; i < nodes.size(); i++) {
        math::Vector3d vi({v[i].X(), v[i].Y(), v[i].Z()});
        math::Point p = ASingLoc + vi;
        p = pl.project(p);
        v[i] = p-ASingLoc;
        v[i].normalize();
    }
    //=====================================================================
    // We get the four vector defining each cross in each node
    std::vector<std::vector<math::Vector3d> > cross_vec;
    cross_vec.resize(nodes.size());
    for (auto i = 0; i < nodes.size(); i++) {
        cross_vec[i].resize(4);
        cross_vec[i][0] = v[i];
        math::AxisAngleRotation rot_pi2(normal, math::Constants::PIDIV2);
        for (auto j = 1; j < 4; j++) {
            cross_vec[i][j] = rot_pi2 * cross_vec[i][j - 1];
        }
    }

    for (auto i = 0; i < nodes.size(); i++) {
        std::vector<double> param_ij;
        int j = (i + 1) % nodes.size();
        Node ni = nodes[i];
        Node nj = nodes[j];
        math::Point pi = ni.point();
        math::Point pj = nj.point();
        double pix = pi.X();
        double piy = pi.Y();
        double piz = pi.Z();
        double pjx = pj.X();
        double pjy = pj.Y();
        double pjz = pj.Z();

        std::vector<math::Vector3d> vis = cross_vec[i];
        std::vector<math::Vector3d> vjs = cross_vec[j];

        bool found = false;
        // We check with the 4 couple of vectors between i and j
        for (auto k = 0; k < 4 && !found; k++) {
            math::Vector3d vi = vis[k];
            math::Vector3d vj;
            double d = -2;
            for (auto l = 0; l < 4; l++) {
                double dot_il = vi.dot(vjs[l]);
                if (dot_il > d) {
                    vj = vjs[l];
                    d = dot_il;
                }
            }

            double vix = vi.X();
            double viy = vi.Y();
            double viz = vi.Z();
            double vjx = vj.X();
            double vjy = vj.Y();
            double vjz = vj.Z();
            // 2nd order equation in alpha the paramrter
            double a = -((nz * piy - ny * piz - nz * pjy + ny * pjz) * vix -
                         (nz * pix - nx * piz - nz * pjx + nx * pjz) * viy +
                         (ny * pix - nx * piy - ny * pjx + nx * pjy) * viz -
                         (nz * piy - ny * piz - nz * pjy + ny * pjz) * vjx +
                         (nz * pix - nx * piz - nz * pjx + nx * pjz) * vjy -
                         (ny * pix - nx * piy - ny * pjx + nx * pjy) * vjz);
            double c = -(cz * ny - cy * nz + nz * piy - ny * piz) * vix +
                       (cz * nx - cx * nz + nz * pix - nx * piz) * viy -
                       (cy * nx - cx * ny + ny * pix - nx * piy) * viz;

            double b = (cz * ny - cy * nz + 2 * nz * piy - 2 * ny * piz - nz * pjy + ny * pjz) * vix -
                       (cz * nx - cx * nz + 2 * nz * pix - 2 * nx * piz - nz * pjx + nx * pjz) * viy +
                       (cy * nx - cx * ny + 2 * ny * pix - 2 * nx * piy - ny * pjx + nx * pjy) * viz -
                       (cz * ny - cy * nz + nz * piy - ny * piz) * vjx +
                       (cz * nx - cx * nz + nz * pix - nx * piz) * vjy -
                       (cy * nx - cx * ny + ny * pix - nx * piy) * vjz;
            std::vector<double> x;
            int nb_x = math::solve2ndDegreePolynomial(a, b, c, x);
            // we can have 0,1 or 2 solutions
            for (auto sol : x) {
                if (sol >= 0 && sol < 1) {
                    param_ij.push_back(sol);
                }
            }  // for(auto sol:x)
        }

        for (auto i = 0; i < param_ij.size(); i++) {
            bool close = false;
            for (auto j = i + 1; j < param_ij.size() && !close; j++) {
                if (abs(param_ij[i] - param_ij[j]) < 1e-2)
                    close = true;
            }
            if (!close) {
                double s = param_ij[i];
                math::Point p_sol = (1 - s) * pi + s * pj;
                math::Vector3d v_sol=p_sol-ASingLoc;
                v_sol.normalize();
                sep.push_back(v_sol);
            }
        }
    }

    // vectors of sep are not necessarily well sorted. We do it now
    std::vector<math::Vector3d> final_sep;
    final_sep.resize(sep.size());

    final_sep[0] = sep[0];
    math::Vector3d ref({sep[0].X(), sep[0].Y(), sep[0].Z()});
    std::vector<double> angle(sep.size(), 0);
    angle[0] = 0;
    for (auto i = 1; i < sep.size(); i++) {
        math::Vector3d vi({sep[i].X(), sep[i].Y(), sep[i].Z()});
        angle[i] = ref.angleIn02PI(vi, normal);
        std::cout << "angle " << i << ": " << angle[i] << std::endl;
    }
    for (auto i = 1; i < sep.size(); i++) {
        int index = 1;
        for (auto j = 1; j < sep.size(); j++) {
            if (i == j)
                continue;
            if (angle[i] > angle[j])
                index++;
        }
        final_sep[index] = sep[i];
        std::cout << index << ": " << final_sep[index] << std::endl;
    }

    return final_sep;
}
