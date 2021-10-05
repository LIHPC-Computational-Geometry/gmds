//
// Created by simon on 12/07/2021.
//

#include "gmds/dualBlocking/EdgeDiscrAlgo.h"
using namespace gmds;
using namespace db;

EdgeDiscrAlgo::EdgeDiscrAlgo(Mesh* AMesh):m_blocks(AMesh) {

}
/*----------------------------------------------------------------------------*/

void EdgeDiscrAlgo::intervalAssignment() {

    for(auto e : m_blocks->edges()){
        Edge edge = m_blocks->get<Edge>(e);

        double length = edge.length();
        int nbSub = ceil(length/0.1);
        targetDiscretization[e] = nbSub;

        double we = 1./length;
        edgeWeight[e] = we;
    }
}
/*----------------------------------------------------------------------------*/
void EdgeDiscrAlgo::init() {

    edgesGrouping();
    intervalAssignment();
    boundaryDiscretization();
}
/*----------------------------------------------------------------------------*/

void EdgeDiscrAlgo::edgesGroupingBlock() {
    //Methode brutale pour obtenir les 3 groupes d'aretes opposées dans chaque bloc
    for(auto b : m_blocks->regions()){
        Region block = m_blocks->get<Region>(b);
        std::vector<Edge> edges = block.get<Edge>();
        std::vector<TCellID> nodes = block.getIDs<Node>();

        //On crée le premier groupe
        //On utilse l'ordre des sommets du bloc, qui eux sont ordonés
        // On va faire 0-1, 2-3, 5-6, 4-7

        TCellID e0,e1,e2,e3;

        for(auto e : edges){
            if((e.getIDs<Node>()[0] == nodes[0] && e.getIDs<Node>()[1] == nodes[1]) ||
                (e.getIDs<Node>()[1] == nodes[0] && e.getIDs<Node>()[0] == nodes[1])){
                e0 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[2] && e.getIDs<Node>()[1] == nodes[3]) ||
                     (e.getIDs<Node>()[1] == nodes[2] && e.getIDs<Node>()[0] == nodes[3])){
                e1 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[4] && e.getIDs<Node>()[1] == nodes[5]) ||
                     (e.getIDs<Node>()[1] == nodes[4] && e.getIDs<Node>()[0] == nodes[6])){
                e2 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[6] && e.getIDs<Node>()[1] == nodes[7]) ||
                     (e.getIDs<Node>()[1] == nodes[6] && e.getIDs<Node>()[0] == nodes[7])){
                e3 = e.id();
            }
        }

        std::tuple<TCellID,TCellID,TCellID,TCellID> groupe1 = std::make_tuple(e0,e1,e2,e3);

        //On crée le premier groupe
        //On utilse l'ordre des sommets du bloc, qui eux sont ordonés
        // On va faire 1-2, 0-3, 4-5, 6-7

        TCellID e4,e5,e6,e7;

        for(auto e : edges){
            if((e.getIDs<Node>()[0] == nodes[1] && e.getIDs<Node>()[1] == nodes[2]) ||
               (e.getIDs<Node>()[1] == nodes[1] && e.getIDs<Node>()[0] == nodes[2])){
                e4 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[0] && e.getIDs<Node>()[1] == nodes[3]) ||
                     (e.getIDs<Node>()[1] == nodes[0] && e.getIDs<Node>()[0] == nodes[3])){
                e5 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[5] && e.getIDs<Node>()[1] == nodes[6]) ||
                     (e.getIDs<Node>()[1] == nodes[5] && e.getIDs<Node>()[0] == nodes[6])){
                e6 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[4] && e.getIDs<Node>()[1] == nodes[7]) ||
                     (e.getIDs<Node>()[1] == nodes[4] && e.getIDs<Node>()[0] == nodes[7])){
                e7 = e.id();
            }
        }

        std::tuple<TCellID,TCellID,TCellID,TCellID> groupe2 = std::make_tuple(e4,e5,e6,e7);

        //On crée le premier groupe
        //On utilse l'ordre des sommets du bloc, qui eux sont ordonés
        // On va faire 1-4, 2-5, 3-6, 0-7

        TCellID e8,e9,e10,e11;

        for(auto e : edges){
            if((e.getIDs<Node>()[0] == nodes[0] && e.getIDs<Node>()[1] == nodes[4]) ||
               (e.getIDs<Node>()[1] == nodes[0] && e.getIDs<Node>()[0] == nodes[4])){
                e8 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[1] && e.getIDs<Node>()[1] == nodes[5]) ||
                     (e.getIDs<Node>()[1] == nodes[1] && e.getIDs<Node>()[0] == nodes[5])){
                e9 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[2] && e.getIDs<Node>()[1] == nodes[6]) ||
                     (e.getIDs<Node>()[1] == nodes[2] && e.getIDs<Node>()[0] == nodes[6])){
                e10 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[3] && e.getIDs<Node>()[1] == nodes[7]) ||
                     (e.getIDs<Node>()[1] == nodes[3] && e.getIDs<Node>()[0] == nodes[7])){
                e11 = e.id();
            }
        }

        std::tuple<TCellID,TCellID,TCellID,TCellID> groupe3 = std::make_tuple(e8,e9,e10,e11);

        std::tuple<std::tuple<TCellID,TCellID,TCellID,TCellID>,std::tuple<TCellID,TCellID,TCellID,TCellID>,std::tuple<TCellID,TCellID,TCellID,TCellID>> block_group = std::make_tuple(groupe1,groupe2,groupe3);

        edgeGroupBlock[b] = block_group;
    }
}
/*----------------------------------------------------------------------------*/

void EdgeDiscrAlgo::boundaryDiscretization() {

    glp_prob *lp;
    int *ia,*ja;
    double *ar;

    int nbCurves = m_blocks->getNbEdges();
    int nbBlocks = m_blocks->getNbRegions();
    int nbFaces = m_blocks->getNbFaces();

    lp = glp_create_prob();
    glp_set_obj_dir(lp,GLP_MIN);
    int nbRows = 2 * nbFaces + nbCurves;
    int nbCols = 3*nbCurves;

    glp_add_rows(lp,nbRows);
    glp_add_cols(lp,nbCols);
    int irow = 1;
    int icol = 1;


    //Constraints
    //this is for face edges equality
    for(int i = 0; i < 2*nbFaces; i++, irow++){
        glp_set_row_bnds(lp,irow,GLP_FX,0,0);
    }
    //this is for De - de - ne = -Ne
    for(int i = 0; i<nbCurves; i++, irow++){
        glp_set_row_name(lp,irow,"D-d-n=-N");
        glp_set_row_bnds(lp,irow,GLP_FX,-targetDiscretization[i],0);
    }


    //Variables
    //this i to enforce a minimum target size
    for(int i = 0;i<nbCurves;i++,icol++){
        std::string ne_name = "ne_"+std::to_string(i);
        glp_set_col_name(lp,icol,ne_name.c_str());
        glp_set_col_kind(lp,icol,GLP_IV);
        glp_set_col_bnds(lp,icol,GLP_LO,1,0);
        glp_set_obj_coef(lp,icol,0);
    }

    //De >= 0
    for(int i = 0;i<nbCurves;i++,icol++){
        std::string De_name = "De_"+std::to_string(i);
        glp_set_col_name(lp,icol,De_name.c_str());
        glp_set_col_kind(lp,icol,GLP_IV);
        glp_set_col_bnds(lp,icol,GLP_LO,0,0);
        glp_set_obj_coef(lp,icol,0);
        glp_set_obj_coef(lp,icol,edgeWeight[i]);
    }

    //de >= 0
    for(int i = 0;i<nbCurves;i++,icol++){
        std::string de_name = "de_"+std::to_string(i);
        glp_set_col_name(lp,icol,de_name.c_str());
        glp_set_col_kind(lp,icol,GLP_IV);
        glp_set_col_bnds(lp,icol,GLP_LO,0,0);
        glp_set_obj_coef(lp,icol,0);
        glp_set_obj_coef(lp,icol,edgeWeight[i]);
    }

    ia = (int*) calloc(1+nbRows*nbCols, sizeof(int));
    ja = (int*) calloc(1+nbRows*nbCols, sizeof(int));
    ar = (double*) calloc(1+nbRows*nbCols, sizeof(double));

    int count = 1;
    for(int i = 0; i<nbRows; i++){
        for(int j = 0; j<nbCols; j++){
            ia[count] = i + 1;
            ja[count] = j + 1;
            count++;
        }
    }

    //fill ar with zeros
    for(int i = 1; i < 1+nbRows*nbCols; i++){
        ar[i] = 0.;
    }


    //fill ar for the block edges equalities
    int indexRow = 1;
    for(int i = 0; i<nbFaces; i++){

        std::tuple<TCellID,TCellID> groupe1 = std::get<0>(edgeGroupFace[i]);
        ar[indexRow + std::get<0>(groupe1)] = 1;
        ar[indexRow + std::get<1>(groupe1)] = -1;
        //ar[indexRow + std::get<2>(groupe1)] = -1;
        //ar[indexRow + std::get<3>(groupe1)] = -1;

        indexRow+=nbCols;

        std::tuple<TCellID,TCellID> groupe2 = std::get<1>(edgeGroupFace[i]);
        ar[indexRow + std::get<0>(groupe2)] = 1;
        ar[indexRow + std::get<1>(groupe2)] = -1;
        //ar[indexRow + std::get<2>(groupe2)] = -1;
        //ar[indexRow + std::get<3>(groupe2)] = -1;

        indexRow+=nbCols;

        /*std::tuple<TCellID,TCellID,TCellID,TCellID> groupe3 = std::get<2>(edgeGroup[i]);
        ar[indexRow + std::get<0>(groupe3)] = 1;
        ar[indexRow + std::get<1>(groupe3)] = 1;
        ar[indexRow + std::get<2>(groupe3)] = -1;
        ar[indexRow + std::get<3>(groupe3)] = -1;

        indexRow+=nbCols;*/
    }
    // De - de - ne =-Ne
    for(int i = 0; i<nbCurves;i++){
        ar[indexRow + i] = -1;
        ar[indexRow + nbCurves + i] = 1;
        ar[indexRow + 2*nbCurves + i] = -1;

        indexRow += nbCols;

    }

    glp_load_matrix(lp,nbRows*nbCols,ia,ja,ar);

    glp_term_out(GLP_ON);

    glp_iocp glpParams;
    glp_init_iocp(&glpParams);
    glpParams.presolve = GLP_ON;

    glp_write_lp(lp,NULL,"block_mapping.txt");

    int glpErr = 0;

    glpErr = glp_intopt(lp,&glpParams);
    switch(glpErr){
        case 0:
            std::cout<<"GLP OK"<<std::endl;
            break;
        case 1:
            std::cout<<"pb solving in GLP"<<std::endl;
            break;
    }

    glpErr = glp_mip_status(lp);
    switch (glpErr){
        case GLP_UNDEF:
            std::cout<<"MIP solution is undefined"<<std::endl;
            throw GMDSException("Block mapping MIP solution is undefined");
            break;
        case GLP_OPT:
            std::cout<<"MIP solution is integer optimal"<<std::endl;
            break;
        case GLP_FEAS:
            std::cout<<"MIP solution is integer feasible"<<std::endl;
            break;
        case GLP_NOFEAS:
            std::cout<<"problem has no integer feasible solution"<<std::endl;
            throw GMDSException("Block mapping problem has no integer feasible solution");
            break;
        default:
            throw GMDSException("Block mapping glp_intopt unknown return code");
    }

    glp_print_mip(lp, "block_mapping_MIP.txt");

    for(int i = 0; i<nbCurves; i++){
        computedDiscretization[i] = glp_mip_col_val(lp,i+1);
    }

    free(ia);
    free(ja);
    free(ar);

    glp_delete_prob(lp);
}
/*----------------------------------------------------------------------------*/

std::map<TCellID, int> EdgeDiscrAlgo::getEdgeDiscretization() {
    return computedDiscretization;
}
/*----------------------------------------------------------------------------*/

void EdgeDiscrAlgo::edgesGrouping() {
    for(auto f : m_blocks->faces()){
        Face face = m_blocks->get<Face>(f);
        std::vector<Edge> edges = face.get<Edge>();
        std::vector<TCellID> nodes = face.getIDs<Node>();

        //On crée le premier groupe
        //On utilse l'ordre des sommets du bloc, qui eux sont ordonés
        // On va faire 0-1, 2-3, 5-6, 4-7

        TCellID e0,e1;

        for(auto e : edges){
            if((e.getIDs<Node>()[0] == nodes[0] && e.getIDs<Node>()[1] == nodes[1]) ||
               (e.getIDs<Node>()[1] == nodes[0] && e.getIDs<Node>()[0] == nodes[1])){
                e0 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[2] && e.getIDs<Node>()[1] == nodes[3]) ||
                     (e.getIDs<Node>()[1] == nodes[2] && e.getIDs<Node>()[0] == nodes[3])){
                e1 = e.id();
            }
        }

        std::tuple<TCellID,TCellID> groupe1 = std::make_tuple(e0,e1);

        //On crée le premier groupe
        //On utilse l'ordre des sommets du bloc, qui eux sont ordonés
        // On va faire 1-2, 0-3

        TCellID e2,e3;

        for(auto e : edges){
            if((e.getIDs<Node>()[0] == nodes[1] && e.getIDs<Node>()[1] == nodes[2]) ||
               (e.getIDs<Node>()[1] == nodes[1] && e.getIDs<Node>()[0] == nodes[2])){
                e2 = e.id();
            }else if((e.getIDs<Node>()[0] == nodes[0] && e.getIDs<Node>()[1] == nodes[3]) ||
                     (e.getIDs<Node>()[1] == nodes[0] && e.getIDs<Node>()[0] == nodes[3])){
                e3 = e.id();
            }
        }

        std::tuple<TCellID,TCellID> groupe2 = std::make_tuple(e2,e3);

        std::tuple<std::tuple<gmds::TCellID,gmds::TCellID>,std::tuple<gmds::TCellID,gmds::TCellID>> edges_groups = std::make_tuple(groupe1,groupe2);

        edgeGroupFace[f] = edges_groups;
    }
}

