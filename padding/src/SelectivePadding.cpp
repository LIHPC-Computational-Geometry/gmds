/*---------------------------------------------------------------------------*/
#include <iostream>
#include <glpk.h>
/*---------------------------------------------------------------------------*/
#include <gmds/padding/SelectivePadding.h>
#include <map>

using namespace gmds;
/*---------------------------------------------------------------------------*/
SelectivePadding::SelectivePadding(Mesh* AMesh)
        : m_mesh(AMesh),
          m_hard_constraint(nullptr),
          m_padding(nullptr),
          m_cplex_with_output(true),
          m_cplex_file_name("prb.lp")
{}
/*---------------------------------------------------------------------------*/
void SelectivePadding::setHardFaces(gmds::Variable<int> *AVar) {
    m_hard_constraint=AVar;
}
/*---------------------------------------------------------------------------*/
void SelectivePadding::setPaddingFaces(gmds::Variable<int> *AVar) {
    m_padding=AVar;
}
/*---------------------------------------------------------------------------*/
void SelectivePadding::setCplexOutputFlag(const bool AWithOutput) {
    m_cplex_with_output=AWithOutput;
}
/*---------------------------------------------------------------------------*/
void SelectivePadding::setCplexFileName(const std::string &AFileName) {
    m_cplex_file_name=AFileName;
}
/*---------------------------------------------------------------------------*/
void SelectivePadding::execute(const Option AOption)
{
    std::cout << "== Execute SelectivePadding =="<< std::endl;
    std::cout<<"Nb cells: "<<m_mesh->getNbRegions()<<std::endl;
    std::cout<<"Nb faces: "<<m_mesh->getNbFaces()<<std::endl;
    std::cout<<"Nb edges: "<<m_mesh->getNbEdges()<<std::endl;
    std::cout<<"Nb nodes: "<<m_mesh->getNbNodes()<<std::endl;

    if(AOption==Option::SBP){


        for(auto edge_id:m_mesh->edges()){
            Edge e=m_mesh->get<Edge>(edge_id);
            std::vector<TCellID> regions_of_e = e.getIDs<Region>();
            if(regions_of_e.size()==1){
                E1H.push_back(edge_id);
            }
            else if(regions_of_e.size()==2){
                E2H.push_back(edge_id);
            }
            else if(regions_of_e.size()==3){
                E3H.push_back(edge_id);
            }
            else if(regions_of_e.size()==4){
                E4H.push_back(edge_id);
            }
        }
        std::cout<<"|E1H|="<<E1H.size()<<std::endl;
        std::cout<<"|E2H|="<<E2H.size()<<std::endl;
        std::cout<<"|E3H|="<<E3H.size()<<std::endl;
        std::cout<<"|E4H|="<<E4H.size()<<std::endl;

        // CALL TO THE SOLVER
        buildAndSolveSBP_pasE1H();
    }
    else{
        std::cout<<"SelectivePadding - Undefined option, so nothing done"
                 <<std::endl;
    }
}
/*---------------------------------------------------------------------------*/
void SelectivePadding::buildAndSolveSBP() {

    std::map<TCellID, int> F2Xi; //faces to Xi index
    std::map<TCellID, int> E2Xi; //edges to Xi index
    //We map integer variables from the solving problem to faces and edges
    // 0 to |edges|-1 are edge ids, then face ids
    int i_x=1;
    for(auto ei_id:m_mesh->edges()){
        E2Xi[ei_id]=i_x;
        i_x++;
    }
    for(auto fi_id:m_mesh->faces()){
        F2Xi[fi_id]=i_x;
        i_x++;
    }

    std::cout<<"E2Xi size "<<E2Xi.size()<<std::endl;
    std::cout<<"F2Xi size "<<F2Xi.size()<<std::endl;


    //==============================================================================
    // Create the problem
    //==============================================================================
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp,"Simple Binary Problem");
    glp_set_obj_dir(lp, GLP_MIN);

    //==============================================================================
    //nb rows = nb auxiliary variables = 1 per edge
    glp_add_rows(lp, m_mesh->getNbEdges());

    for (auto ei_id:m_mesh->edges()) {
        std::string aux_name = std::string("E")+std::to_string(ei_id);
        glp_set_row_name(lp, E2Xi[ei_id] , aux_name.c_str()); //name
        glp_set_row_bnds(lp, E2Xi[ei_id] , GLP_FX,0,0);          //Binary variable
    }

    //==============================================================================
    // nb columns = nb structural variables = 1 per face + 1 per edge
    double coeff = pow((double)(m_mesh->getNbHexahedra()), (-2.0 / 3.0));

    //std::cout<<"|H|^{-2/3}|="<<coeff<<std::endl;

    glp_add_cols(lp, m_mesh->getNbEdges() + m_mesh->getNbFaces());

    for (auto ei_id:m_mesh->edges()) {
        std::string name = "e"+std::to_string(ei_id);
        glp_set_col_name(lp,E2Xi[ei_id] , name.c_str());
        glp_set_col_kind(lp, E2Xi[ei_id] , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi[ei_id] , 0); //no impact on Energy value
    }

    // Epadding = |H|^{-2/3} sum Xfi

    for (auto fi_id:m_mesh->faces()) {
        if (m_hard_constraint->value(fi_id) == 1) {
            //HARD CONSTRAINT = 1
            std::string name = "HF"+std::to_string(fi_id);
            glp_set_col_name(lp, F2Xi[fi_id] , name.c_str());
            glp_set_col_bnds(lp, F2Xi[fi_id] , GLP_FX, 1.0, 1.0);
            glp_set_obj_coef(lp, F2Xi[fi_id] , coeff);
        }
        else {
            std::string name = "f"+std::to_string(fi_id);
            glp_set_col_name(lp, F2Xi[fi_id] , name.c_str());
            glp_set_col_kind(lp, F2Xi[fi_id] , GLP_BV);          //Binary variable
            glp_set_obj_coef(lp, F2Xi[fi_id] , coeff);
        }
    }

    //==============================================================================
    // Now we build the matrix of constraints w
    //TODO check those values that are totally empirical right now
    int    ia[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];
    int    ja[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];
    double ar[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];

    //ia row number, so an edge id
    //ja column number, so a face id
    //ar weight associated to the face_id always 1/2 for us
    //Now we build the matrix of constraints
    int i=1;
    for (auto ei_id:m_mesh->edges()) {
        Edge ei = m_mesh->get<Edge>(ei_id);
        ia[i] = E2Xi[ei_id];
        ja[i] = E2Xi[ei_id];
        ar[i]= -2;
        i++;
        std::vector<TCellID > adj_face_ids = ei.getIDs<Face>();
        for(auto fi_id:adj_face_ids){
            ia[i] =E2Xi[ei_id];
            ja[i] =F2Xi[fi_id];
            ar[i]= 1;
            i++;
        }
    }

    glp_load_matrix(lp, i-1 , ia, ja, ar);

    if(m_cplex_with_output)   {
        glp_write_lp(lp,0,m_cplex_file_name.c_str());
    }
    //std::cout<<"Nb binary var= "<<glp_get_num_bin(lp)<<std::endl;

    int err= glp_simplex(lp,NULL);
    glp_intopt(lp,0);//&param);

    double E_Padding = glp_mip_obj_val(lp);

  //  std::string sol_file = "lp_sol";
  //  glp_write_sol(lp,sol_file.c_str());

    //Getting the function of extraction des valeurs de Xfi pour les affichées dans paraview
    int HF=0, F_HF=0;
    for(auto fi_id:m_mesh->faces()){
        double val = glp_mip_col_val(lp, F2Xi[fi_id] );
        (*m_padding)[fi_id] = val;
        if(val==0)
            F_HF++;
        else
            HF++;
    }

    int e=0;
    for(auto ei_id:m_mesh->edges()){
        double val = glp_mip_col_val(lp, E2Xi[ei_id] );;
        if(val!=0)
            e++;
    }

    std::cout << "Energy ="<<E_Padding<<std::endl;
    std::cout << "|HF|   =" << HF << std::endl;
    std::cout << "|F-HF| =" << F_HF << std::endl;
    std::cout << "|E|    =" << e << std::endl;
    glp_delete_prob(lp);
}
/*----------------------------------------------------------------*/
void SelectivePadding::buildAndSolveSBP_pasE1H() {

    std::vector<TCellID> HF;
    for(auto f:m_mesh->faces()){
        if(m_hard_constraint->value(f) == 1){
            HF.push_back(f);
        }
    }

    std::map<TCellID, int> F2Xi; //faces to Xi index
    std::map<TCellID, int> E2Xi; //edges to Xi index
    std::map<TCellID, int> HF2Xi;

    std::map<TCellID,std::pair<TCellID,TCellID >> tE2H;

    // tE3H_E4H[eid] -> <(fid1,fid2),(fid3,fid4)> eid is an edge of valence 3 or 4 and the value is a pair of pairs of
    // parallel faces
    std::map<TCellID,std::pair<TCellID,TCellID >> tE3H_4H;

    for(auto e : E2H){
        Edge edge = m_mesh->get<Edge>(e);
        std::vector<Face> e_faces = edge.get<Face>();
        if(e_faces[0].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[1].id(),e_faces[2].id());
        }else if(e_faces[1].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[0].id(),e_faces[2].id());
        }else if(e_faces[2].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[0].id(),e_faces[1].id());
        }
    }
    std::cout<<"tE2H "<<tE2H.size()<<std::endl;
    for(auto e : E3H){
        Edge edge = m_mesh->get<Edge>(e);
        std::pair<TCellID,TCellID> pair1;
        bool first = true;
        std::pair<TCellID,TCellID> pair2;
        std::vector<Face> e_faces = edge.get<Face>();
        for(auto f : e_faces){
            if(f.get<Region>().size() == 1){
                int r_id = f.getIDs<Region>()[0];
                for(auto f2 : e_faces){
                    if(f2.get<Region>().size() == 2){
                        if(f2.getIDs<Region>()[0] != r_id && f2.getIDs<Region>()[1] != r_id){
                            if(first) {
                                pair1 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                //first = false;
                                break;
                            }
                            else{

                                pair2 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                break;
                            }
                        }
                    }
                }
            }
        }
        tE3H_4H[e] = pair1;
    }

    for(auto e : E4H){
        Edge edge = m_mesh->get<Edge>(e);
        std::vector<TCellID> faces_used;
        std::pair<TCellID,TCellID> pair1;
        int first = 0;
        std::pair<TCellID,TCellID> pair2;
        std::vector<Face> e_faces = edge.get<Face>();
        for(auto f : e_faces){
            if(std::find(faces_used.begin(),faces_used.end(),f.id()) == faces_used.end()){
                faces_used.push_back(f.id());
                int r_id0 = f.getIDs<Region>()[0];
                int r_id1 = f.getIDs<Region>()[1];
                for(auto f2 : e_faces){
                    if(f2.get<Region>().size() == 2){
                        if((f2.getIDs<Region>()[0] != r_id0 && f2.getIDs<Region>()[1] != r_id1) &&
                           (f2.getIDs<Region>()[1] != r_id0 && f2.getIDs<Region>()[0] != r_id1)){
                            if(first == 0) {
                                pair1 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                first = 1;
                                faces_used.push_back(f2.id());
                                break;
                            }
                            else if(first == 1){
                                pair2 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                faces_used.push_back(f2.id());
                                first = 2;
                            }
                        }
                    }
                }
                if(first == 1) break;
            }
        }
        tE3H_4H[e] = pair1;
    }
    //We map integer variables from the solving problem to faces and edges
    // 0 to |edges|-1 are edge ids, then face ids


    int i_x=1;

    for(auto ei_id:E3H){
        E2Xi[ei_id]=i_x;
        i_x++;
    }
    for(auto ei_id:E4H){
        E2Xi[ei_id]=i_x;
        i_x++;
    }

    i_x=1;
    std::map<int,std::pair<TCellID,TCellID>> tF;
    for(auto pair : tE2H){
        tF[i_x] = pair.second;
        i_x++;
    }

    for(auto pair : tE3H_4H){
        tF[i_x] = pair.second;
        i_x++;
    }

    i_x=1;
    for(auto fi_id:m_mesh->faces()){
        F2Xi[fi_id]=i_x;
        i_x++;
    }

    std::cout<<"E2Xi size "<<E2Xi.size()<<std::endl;
    std::cout<<"F2Xi size "<<F2Xi.size()<<std::endl;


    //==============================================================================
    // Create the problem
    //==============================================================================
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp,"Simple Binary Problem");
    glp_set_obj_dir(lp, GLP_MIN);

    //==============================================================================
    //nb rows = nb auxiliary variables = 1 per edge
    glp_add_rows(lp, tF.size() + E2Xi.size());


    for(auto ei:E3H){
        std::string aux_name = std::string("E")+std::to_string(ei);
        glp_set_row_name(lp, E2Xi[ei], aux_name.c_str()); //name
        glp_set_row_bnds(lp, E2Xi[ei], GLP_FX,0,0);          //Binary variable
    }
    for (auto pair:tF) {
        std::string aux_name = std::string("tF")+std::to_string(pair.first);
        glp_set_row_name(lp, E2Xi.size() + pair.first, aux_name.c_str()); //name
        glp_set_row_bnds(lp, E2Xi.size() + pair.first, GLP_FX,0,0);          //Binary variable
    }

    //==============================================================================
    // nb columns = nb structural variables = 1 per face + 1 per edge
    double coeff_Padd = pow((double)(m_mesh->getNbHexahedra()), (-2.0 / 3.0));

    double coeff_Comp = pow((double)(m_mesh->getNbHexahedra()), (-1.0 / 3.0));

    //std::cout<<"|H|^{-2/3}|="<<coeff<<std::endl;

    glp_add_cols(lp, E2Xi.size() + (2*tF.size()) + m_mesh->getNbFaces());

    for(auto ei:E3H){
        std::string aux_name = std::string("E")+std::to_string(ei);
        glp_set_col_name(lp, E2Xi[ei], aux_name.c_str()); //name
        glp_set_col_kind(lp, E2Xi[ei], GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi[ei], 0);
    }
    for(auto ei:E4H){
        std::string aux_name = std::string("E")+std::to_string(ei);
        glp_set_col_name(lp, E2Xi[ei], aux_name.c_str()); //name
        glp_set_col_kind(lp, E2Xi[ei], GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi[ei], 0);
    }
    for (auto pair:tF) {
        std::string name = "Xtf"+std::to_string(pair.first)+std::string("_1");
        glp_set_col_name(lp, E2Xi.size() + pair.first , name.c_str());
        glp_set_col_kind(lp, E2Xi.size() + pair.first , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi.size() + pair.first , coeff_Comp);
    }
    for (auto pair:tF) {
        std::string name = "Xtf"+std::to_string(pair.first)+std::string("_2");
        glp_set_col_name(lp, E2Xi.size() + tF.size() + pair.first , name.c_str());
        glp_set_col_kind(lp, E2Xi.size() + tF.size() + pair.first , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi.size() + tF.size() + pair.first , coeff_Comp);
    }

    // Epadding = |H|^{-2/3} sum Xfi

    for (auto fi_id:m_mesh->faces()) {
        if (m_hard_constraint->value(fi_id) == 1) {
            //HARD CONSTRAINT = 1
            std::string name = "HF"+std::to_string(fi_id);
            glp_set_col_name(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , name.c_str());
            glp_set_col_bnds(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , GLP_FX, 1.0, 1.0);
            glp_set_obj_coef(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , coeff_Padd);
        }
        else {
            std::string name = "f"+std::to_string(fi_id);
            glp_set_col_name(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , name.c_str());
            glp_set_col_kind(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , GLP_BV);          //Binary variable
            glp_set_obj_coef(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] , coeff_Padd);
        }
    }

    //==============================================================================
    // Now we build the matrix of constraints w
    //TODO check those values that are totally empirical right now
    int    ia[1 + 4*m_mesh->getNbEdges() + 100*m_mesh->getNbFaces()];
    int    ja[1 + 4*m_mesh->getNbEdges() + 100*m_mesh->getNbFaces()];
    double ar[1 + 4*m_mesh->getNbEdges() + 100*m_mesh->getNbFaces()];

    //ia row number, so an edge id
    //ja column number, so a face id
    //ar weight associated to the face_id always 1/2 for us
    //Now we build the matrix of constraints
    int i=1;
    for (auto ei_id:E3H) {
        Edge ei = m_mesh->get<Edge>(ei_id);
        ia[i] = E2Xi[ei_id];
        ja[i] = E2Xi[ei_id];
        ar[i]= -2;
        i++;
        std::vector<TCellID > adj_face_ids = ei.getIDs<Face>();
        for(auto fi_id:adj_face_ids){
            ia[i] =E2Xi[ei_id];
            ja[i] =E2Xi.size() + (2*tF.size()) + F2Xi[fi_id];
            ar[i]= 1;
            i++;
        }
    }
    for (auto ei_id:E4H) {
        Edge ei = m_mesh->get<Edge>(ei_id);
        ia[i] = E2Xi[ei_id];
        ja[i] = E2Xi[ei_id];
        ar[i]= -2;
        i++;
        std::vector<TCellID > adj_face_ids = ei.getIDs<Face>();
        for(auto fi_id:adj_face_ids){
            ia[i] =E2Xi[ei_id];
            ja[i] =E2Xi.size() + (2*tF.size()) + F2Xi[fi_id];
            ar[i]= 1;
            i++;
        }
    }
    for (auto pair:tF) {
        ia[i] =E2Xi.size() + pair.first;
        ja[i] = E2Xi.size() + (2*tF.size()) + F2Xi[pair.second.first];
        ar[i] = 1;
        std::cout<<i<<" face1 : "<<ia[i]<<","<<ja[i]<<" = "<<ar[i]<<"\n"<<std::endl;
        std::cout<<"face = "<<pair.second.first<<std::endl;
        /*std::cout<<"ia["<<i<<"] = "<<ia[i]<<std::endl;
        std::cout<<"ja["<<i<<"] = "<<ja[i]<<std::endl;
        std::cout<<"ar["<<i<<"] = "<<ar[i]<<std::endl;
        std::cout<<std::endl;*/
        i++;

        ia[i] = E2Xi.size() + pair.first;
        ja[i] = E2Xi.size() + (2*tF.size()) + F2Xi[pair.second.second];
        ar[i] = -1;
        std::cout<<i<<" face2 : "<<ia[i]<<","<<ja[i]<<" = "<<ar[i]<<"\n"<<std::endl;
        std::cout<<"face = "<<pair.second.second<<std::endl;
        /*std::cout<<"ia["<<i<<"] = "<<ia[i]<<std::endl;
        std::cout<<"ja["<<i<<"] = "<<ja[i]<<std::endl;
        std::cout<<"ar["<<i<<"] = "<<ar[i]<<std::endl;
        std::cout<<std::endl;*/
        i++;

        ia[i] = E2Xi.size() + pair.first;
        ja[i] = E2Xi.size() + pair.first;
        ar[i] = -1;
        std::cout<<i<<" var1 : "<<ia[i]<<","<<ja[i]<<" = "<<ar[i]<<"\n"<<std::endl;
        /*std::cout<<"ia["<<i<<"] = "<<ia[i]<<std::endl;
        std::cout<<"ja["<<i<<"] = "<<ja[i]<<std::endl;
        std::cout<<"ar["<<i<<"] = "<<ar[i]<<std::endl;
        std::cout<<std::endl;*/
        i++;

        ia[i] = E2Xi.size() + pair.first;
        ja[i] = E2Xi.size() + tF.size() + pair.first;
        ar[i] = 1;
        std::cout<<i<<" var2 : "<<ia[i]<<","<<ja[i]<<" = "<<ar[i]<<"\n"<<std::endl;
        /*std::cout<<"ia["<<i<<"] = "<<ia[i]<<std::endl;
        std::cout<<"ja["<<i<<"] = "<<ja[i]<<std::endl;
        std::cout<<"ar["<<i<<"] = "<<ar[i]<<std::endl;
        std::cout<<std::endl;*/
        i++;

    }

    glp_load_matrix(lp, i-1 , ia, ja, ar);

    if(m_cplex_with_output)   {
        glp_write_lp(lp,0,m_cplex_file_name.c_str());
    }
    //std::cout<<"Nb binary var= "<<glp_get_num_bin(lp)<<std::endl;

    int err= glp_simplex(lp,NULL);
    glp_intopt(lp,0);//&param);

    double E_Padding = glp_mip_obj_val(lp);

      std::string sol_file = "lp_sol";
      glp_write_sol(lp,sol_file.c_str());

    //Getting the function of extraction des valeurs de Xfi pour les affichées dans paraview
    int nbHF=0, F_HF=0;
    for(auto fi_id:m_mesh->faces()){
        double val = glp_mip_col_val(lp, E2Xi.size() + (2*tF.size()) + F2Xi[fi_id] );
        (*m_padding)[fi_id] = val;
        if(val==0)
            F_HF++;
        else
            nbHF++;
    }

    int e=0;
    for(auto ei:E3H){
        double val = glp_mip_col_val(lp, E2Xi[ei]);
        if(val!=0)
            e++;
    }
    for(auto ei:E4H){
        double val = glp_mip_col_val(lp, E2Xi[ei]);
        if(val!=0)
            e++;
    }
    for(int i_c = 1; i_c<(2*tF.size()) + m_mesh->getNbFaces();i_c++){
        double val = glp_mip_col_val(lp, i_c);
        //std::cout<<"val ["<<i_c<<"] = "<<val<<std::endl;

    }

    std::cout << "Energy ="<<E_Padding<<std::endl;
    std::cout << "|HF|   =" << nbHF << std::endl;
    std::cout << "|F-HF| =" << F_HF << std::endl;
    std::cout << "|E|    =" << e << std::endl;
    glp_delete_prob(lp);
}
/*----------------------------------------------------------------*/
void SelectivePadding::buildAndSolveSBP_test() {

    std::map<TCellID, int> F2Xi; //faces to Xi index
    std::map<TCellID, int> E2Xi; //edges to Xi index

    std::vector<TCellID> E3H_E4H;
    for(auto e : E3H){
        E3H_E4H.push_back(e);
    }
    for(auto e : E4H){
        E3H_E4H.push_back(e);
    }

    // Map for the turn constraint
    // tE2H[eid] -> (fid1,fid2) eid is an edge of valence 2 and fid1 and fid2 are the faces on the boundary (ie. parallels)
    std::map<TCellID,std::pair<TCellID,TCellID >> tE2H;

    // tE3H_E4H[eid] -> <(fid1,fid2),(fid3,fid4)> eid is an edge of valence 3 or 4 and the value is a pair of pairs of
    // parallel faces
    std::map<TCellID,std::pair<std::pair<TCellID,TCellID >,std::pair<TCellID,TCellID >>> tE3H_E4H;

    for(auto e : E2H){
        Edge edge = m_mesh->get<Edge>(e);
        std::vector<Face> e_faces = edge.get<Face>();
        if(e_faces[0].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[1].id(),e_faces[2].id());
        }else if(e_faces[1].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[0].id(),e_faces[2].id());
        }else if(e_faces[2].getIDs<Region>().size() == 2){
            tE2H[e] = std::pair<TCellID,TCellID>(e_faces[0].id(),e_faces[1].id());
        }
    }
    for(auto e : E3H){
        Edge edge = m_mesh->get<Edge>(e);
        std::pair<TCellID,TCellID> pair1;
        bool first = true;
        std::pair<TCellID,TCellID> pair2;
        std::vector<Face> e_faces = edge.get<Face>();
        for(auto f : e_faces){
            if(f.get<Region>().size() == 1){
                int r_id = f.getIDs<Region>()[0];
                for(auto f2 : e_faces){
                    if(f2.get<Region>().size() == 2){
                        if(f2.getIDs<Region>()[0] != r_id && f2.getIDs<Region>()[1] != r_id){
                            if(first) {
                                pair1 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                first = false;
                            }
                            else{

                                pair2 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                break;
                            }
                        }
                    }
                }
            }
        }
        tE3H_E4H[e] = std::pair<std::pair<TCellID,TCellID >,std::pair<TCellID,TCellID >>(pair1,pair2);
    }
    for(auto e : E4H){
        Edge edge = m_mesh->get<Edge>(e);
        std::vector<TCellID> faces_used;
        std::pair<TCellID,TCellID> pair1;
        int first = 0;
        std::pair<TCellID,TCellID> pair2;
        std::vector<Face> e_faces = edge.get<Face>();
        for(auto f : e_faces){
            if(std::find(faces_used.begin(),faces_used.end(),f.id()) == faces_used.end()){
                faces_used.push_back(f.id());
                int r_id0 = f.getIDs<Region>()[0];
                int r_id1 = f.getIDs<Region>()[1];
                for(auto f2 : e_faces){
                    if(f2.get<Region>().size() == 2){
                        if((f2.getIDs<Region>()[0] != r_id0 && f2.getIDs<Region>()[1] != r_id1) &&
                        (f2.getIDs<Region>()[1] != r_id0 && f2.getIDs<Region>()[0] != r_id1)){
                            if(first == 0) {
                                pair1 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                first = 1;
                                faces_used.push_back(f2.id());
                            }
                            else if(first == 1){
                                pair2 = std::pair<TCellID,TCellID>(f.id(),f2.id());
                                faces_used.push_back(f2.id());
                                first = 2;
                            }
                        }
                    }
                }
                if(first == 2) break;
            }
        }
        tE3H_E4H[e] = std::pair<std::pair<TCellID,TCellID >,std::pair<TCellID,TCellID >>(pair1,pair2);
    }

    Variable<int>* var_pair = m_mesh->newVariable<int, GMDS_FACE>("pair");


    std::cout<<"tE3H4H size "<<tE3H_E4H.size()<<std::endl;
    //ici techniquement on a le bon nombre de paires de faces
    //We map integer variables from the solving problem to faces and edges
    // 0 to |edges|-1 are edge ids, then face ids
    int i_x=1;
    for(auto ei_id:E3H){
        E2Xi[ei_id]=i_x;
        i_x++;
    }
    for(auto ei_id:E4H){
        E2Xi[ei_id]=i_x;
        i_x++;
    }


    int i_tE = i_x; // Pourquoi ? Parce que i_x nous donne le nombre de rows des contraintes sur les arêtes
    // l'indice de départ pour les contraintes sur les paires de faces commence donc au i_x ieme
    std::map<int,std::pair<TCellID,TCellID>> tF;
    for(auto pair : tE2H){
        tF[i_tE] = pair.second;
        i_tE++;
    }

    for(auto pair : tE3H_E4H){
        tF[i_tE] = pair.second.first;
        i_tE++;
        tF[i_tE] = pair.second.second;
        i_tE++;
    }

    for(auto pair : tF){
        if(var_pair->value(pair.second.first) == 0){
            var_pair->set(pair.second.first, pair.first);
        }

        if(var_pair->value(pair.second.second) == 0){
            var_pair->set(pair.second.second, pair.first);
        }
    }

    // Comme on a utilisé un autre indice pour les contraintes de paires, on peut continuer avec i_x pour les colonnes
    // pour les variables des faces
    for(auto fi_id:m_mesh->faces()){
        F2Xi[fi_id]=i_tE;
        i_tE++;
    }

    std::cout<<"E2Xi size "<<E2Xi.size()<<std::endl;
    std::cout<<"F2Xi size "<<F2Xi.size()<<std::endl;



    //==============================================================================
    // Create the problem
    //==============================================================================
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp,"Simple Binary Problem");
    glp_set_obj_dir(lp, GLP_MIN);

    //==============================================================================
    //nb rows = one per E2H edges for the turn constraint, 2 per E3H_E4H for the 2 pair of parallel edges for the turn
    //constraint and one per E3H_E4H edge for the padding constraint
    glp_add_rows(lp, E3H_E4H.size() + tF.size());

    // ça c'est ok
    for (auto ei_id:E3H_E4H) {
        std::string aux_name = std::string("E")+std::to_string(ei_id);
        glp_set_row_name(lp, E2Xi[ei_id] , aux_name.c_str()); //name
        glp_set_row_bnds(lp, E2Xi[ei_id] , GLP_FX,0,0);          //Binary variable
    }


    for (auto pair:tF) {
        //std::string aux_name = std::string("tf")+std::to_string(pair.second.first)+std::string("_f")+std::to_string(pair.second.second);
        std::string aux_name = std::string("TF")+std::to_string(pair.first);
        glp_set_row_name(lp, pair.first , aux_name.c_str()); //name
        glp_set_row_bnds(lp, pair.first , GLP_FX,0,0);          //Binary variable
    }

    //==============================================================================
    // nb columns = nb structural variables = 1 per face + 1 per edge
    double coeff = pow((double)(m_mesh->getNbHexahedra()), (-2.0 / 3.0));

    //std::cout<<"|H|^{-2/3}|="<<coeff<<std::endl;

    glp_add_cols(lp, E3H_E4H.size() + tF.size() + m_mesh->getNbFaces());

    for (auto ei_id:E3H_E4H) {
        std::string name = "e"+std::to_string(ei_id);
        glp_set_col_name(lp, E2Xi[ei_id] , name.c_str());
        glp_set_col_kind(lp, E2Xi[ei_id] , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi[ei_id] , 0); //no impact on Energy value
    }

    for (auto tf:tF) {
        std::string name = "tf"+std::to_string(tf.first);
        glp_set_col_name(lp, tf.first , name.c_str());
        glp_set_col_kind(lp, tf.first , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, tf.first , 0); //no impact on Energy value
    }

    // Epadding = |H|^{-2/3} sum Xfi

    for (auto fi_id:m_mesh->faces()) {
        if (m_hard_constraint->value(fi_id) == 1) {
            //HARD CONSTRAINT = 1
            std::string name = "HF"+std::to_string(fi_id);
            glp_set_col_name(lp, F2Xi[fi_id] , name.c_str());
            glp_set_col_bnds(lp, F2Xi[fi_id] , GLP_FX, 1.0, 1.0);
            glp_set_obj_coef(lp, F2Xi[fi_id] , coeff);
        }
        else {
            std::string name = "f"+std::to_string(fi_id);
            glp_set_col_name(lp, F2Xi[fi_id] , name.c_str());
            glp_set_col_kind(lp, F2Xi[fi_id] , GLP_BV);          //Binary variable
            glp_set_obj_coef(lp, F2Xi[fi_id] , coeff);
        }
    }

    //==============================================================================
    // Now we build the matrix of constraints w
    //TODO check those values that are totally empirical right now
    int    ia[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];
    int    ja[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];
    double ar[1 + 4*m_mesh->getNbEdges() + 4*m_mesh->getNbFaces()];

    //ia row number, it's the number of constraints
    //ja column number, so a face id
    //ar weight associated to the face_id always 1/2 for us
    //Now we build the matrix of constraints
    int i=1;

    for (auto ei_id:E3H_E4H) {
        Edge ei = m_mesh->get<Edge>(ei_id);
        ia[i] = E2Xi[ei_id];
        ja[i] = E2Xi[ei_id];
        ar[i]= -2;
        i++;
        std::vector<TCellID > adj_face_ids = ei.getIDs<Face>();
        for(auto fi_id:adj_face_ids){
            ia[i] =E2Xi[ei_id];
            ja[i] =F2Xi[fi_id];
            ar[i]= 1;
            i++;
        }
    }
    for(auto pair:tF){

        ia[i] = pair.first;
        ja[i] = pair.first;
        ar[i] = -2;
        i++;

        //std::cout<<"ia["<<i<<"] = "<<pair.first<<"; ja["<<i<<"] = "<<F2Xi[pair.second.first]<<std::endl;
        ia[i] = pair.first;
        ja[i] = F2Xi[pair.second.first];
        ar[i]= 1;
        i++;

        //std::cout<<"ia["<<i<<"] = "<<pair.first<<"; ja["<<i<<"] = "<<F2Xi[pair.second.second]<<std::endl;
        ia[i] = pair.first;
        ja[i] = F2Xi[pair.second.second];
        ar[i]= 1;
        i++;
    }

    glp_load_matrix(lp, i-1 , ia, ja, ar);

    if(m_cplex_with_output)   {
        glp_write_lp(lp,0,m_cplex_file_name.c_str());
    }
    //std::cout<<"Nb binary var= "<<glp_get_num_bin(lp)<<std::endl;

    int err= glp_simplex(lp,NULL);
    glp_intopt(lp,0);//&param);

    double E_Padding = glp_mip_obj_val(lp);

    //  std::string sol_file = "lp_sol";
    //  glp_write_sol(lp,sol_file.c_str());

    //Getting the function of extraction des valeurs de Xfi pour les affichées dans paraview
    int HF=0, F_HF=0;
    for(auto fi_id:m_mesh->faces()){
        double val = glp_mip_col_val(lp, F2Xi[fi_id] );
        (*m_padding)[fi_id] = val;
        if(val==0)
            F_HF++;
        else
            HF++;
    }

    int e=0;
    for(auto ei_id:E3H_E4H){
        double val = glp_mip_col_val(lp, E2Xi[ei_id] );;
        if(val!=0)
            e++;
    }

    std::cout << "Energy ="<<E_Padding<<std::endl;
    std::cout << "|HF|   =" << HF << std::endl;
    std::cout << "|F-HF| =" << F_HF << std::endl;
    std::cout << "|E|    =" << e << std::endl;
    glp_delete_prob(lp);
}