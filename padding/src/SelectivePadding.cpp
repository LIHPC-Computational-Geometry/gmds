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
          m_cplex_with_output(false),
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
        std::vector<TCellID> E1H, E2H, E3H;

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
        }
        std::cout<<"|E1H|="<<E1H.size()<<std::endl;
        std::cout<<"|E2H|="<<E2H.size()<<std::endl;
        std::cout<<"|E3H|="<<E3H.size()<<std::endl;

        // CALL TO THE SOLVER
        buildAndSolveSBP();
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

    glp_add_cols(lp, m_mesh->getNbEdges()+m_mesh->getNbFaces());

    for (auto ei_id:m_mesh->edges()) {
        std::string name = "e"+std::to_string(ei_id);
        glp_set_col_name(lp,E2Xi[ei_id] , name.c_str());
        glp_set_col_kind(lp, E2Xi[ei_id] , GLP_BV);          //Binary variable
        glp_set_obj_coef(lp, E2Xi[ei_id] , 0); //no impact on Energy value
    }

    for (auto fi_id:m_mesh->faces()) {
        if (m_hard_constraint->value(fi_id) == 1) {
            //HARD CONSTRAINT = 1
            std::string name = "hf"+std::to_string(fi_id);
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
