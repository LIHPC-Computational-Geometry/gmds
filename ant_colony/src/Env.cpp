/******************************************************************************/
#include <gmds/ant_colony/Env.h>
/******************************************************************************/
using namespace gmds;
using namespace math;

/******************************************************************************/


std::string format_number(int number, int max_digit) {
    std::string s;
    int size, i, diff;
    s = "";
    size = std::to_string(number).size();
    diff = max_digit - size;
    i = 0;
    while (i < diff) {
        s += "0";
        i++;
    }
    return s + std::to_string(number);
}

void Env::writeVTK(std::string fname) {
    if(this->write_vtk) {
        IGMeshIOService ios(&(this->m));
        VTKWriter writer(&ios);
        writer.setCellOptions(N | R);
        writer.setDataOptions(N | R);
        writer.write(fname); // paraview ctrl + space (Shrink) pour avoir uniquement les hex (enlever les F dans l'écriture)
    }
}

void Env::writeSolution(std::vector<char> solution, std::string fname,
                        std::vector<float> pheromones,std::vector<TCellID> history_faces, std::vector<TCellID> history_edges) {


    Variable<int> * history_faces_visu = m.newVariable<int, gmds::GMDS_FACE>("history_faces");
    Variable<int> * history_edges_visu = m.newVariable<int, gmds::GMDS_EDGE>("history_edges");
    history_faces_visu->setValuesTo(0);
    history_edges_visu->setValuesTo(0);
    int cpt = 0;
    for(auto x : history_faces){
        history_faces_visu->set(x,cpt);
        cpt++;
    }
    cpt=0;
    for(auto x : history_edges){
        history_edges_visu->set(x,cpt);
        cpt++;
    }
    for (auto it = solution.begin(); it != solution.end(); ++it) {
        std::cout << std::distance(solution.begin(), it) << ' ' << *it << std::endl;
        this->face_in_solution->set(std::distance(solution.begin(), it), int(*it));
    }
    std::cout << std::endl;

    for (auto x: this->initial_faces) {
        this->face_in_solution->set(x, FaceStatus::Initial);
    }
    writeVTK(fname);
    resetVariables();
    m.deleteVariable(gmds::GMDS_FACE,"history_faces");
    m.deleteVariable(gmds::GMDS_EDGE,"history_edges");
}

Env::Env(TInt x_n, TInt y_n, TInt z_n, const int *exist_flatten) : m(
        gmds::Mesh(gmds::MeshModel(
                DIM3 | R | F | E | N | R2N | R2F | R2E | F2N | F2R | F2E | E2F | E2N | N2E | N2R | E2R))),
                                                                   x_n(x_n), y_n(y_n), z_n(z_n) {
    std::cout << "Env initialized" << std::endl;
    // GridBuilder
    gmds::GridBuilder gb(&m, 3);
    // Number of bloc for each dim
    this->x_n++;
    this->y_n++;
    this->z_n++;
    gb.execute(this->x_n, 1.0, this->y_n, 1.0, this->z_n, 1.0);
    std::cout << m.getNbRegions() << std::endl;
    assert(m.getNbRegions() == (this->x_n - 1) * (this->y_n - 1) * (this->z_n - 1));
    gmds::Region r_tmp;
    // Building our gridlike structure
    // Assign boolean value to blocs (exists or not) 0 no 1 yes
    //this->bloc_exist = this->m.newVariable<int, gmds::GMDS_REGION>("bloc_exist");
    //bloc_exist->setValuesTo(0);
    this->var_delete = this->m.newVariable<int, gmds::GMDS_NODE>("var_delete");
    var_delete->setValuesTo(0);

    for (auto i: m.regions()) {
        if (exist_flatten[i] == 0) {
            m.deleteRegion(i);
        }
    }
    for (auto i: m.regions()) {
        gmds::Region r = m.get<gmds::Region>(i);
        std::vector<TCellID> n_ids = r.getIDs<Node>();
        for (auto n: n_ids) {
            var_delete->set(n, 1);
        }
    }

    std::cout << "Nodes " << m.getNbNodes() << std::endl;
    std::cout << "Regions " << m.getNbRegions() << std::endl;
    std::cout << "Edges " << m.getNbEdges() << std::endl;
    std::cout << "Faces " << m.getNbFaces() << std::endl;
    std::cout << "Hexaedra " << m.getNbHexahedra() << std::endl;

    // Correct numberof bloc
    this->x_n--, this->y_n--, this->z_n--;
    std::cout << "(X, Y, Z) = ( " << x_n << ", " << y_n << ", " << z_n << ")" << std::endl;
    //BoundaryOperator bndo(&m);
    // Mesh doctor to generate connectivity
    gmds::MeshDoctor doc(&m);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    std::cout << "Nodes " << m.getNbNodes() << std::endl;
    std::cout << "Regions " << m.getNbRegions() << std::endl;
    std::cout << "Edges " << m.getNbEdges() << std::endl;
    std::cout << "Faces " << m.getNbFaces() << std::endl;
    std::cout << "Hexaedra " << m.getNbHexahedra() << std::endl;

    this->edge_boundary = m.newVariable<int, gmds::GMDS_EDGE>("edge_boundary");
    this->face_boundary = m.newVariable<int, gmds::GMDS_FACE>("face_boundary");

    // Faces au bord
    for (auto f: m.faces()) {
        gmds::Face mface = m.get<gmds::Face>(f);
        std::vector<TCellID> v_region;
        mface.delegateGetRegionIDs(v_region);
        if (v_region.size() == 1) {
            face_boundary->set(f, 1);
        }
    }
    // Arêtes au bord
    for (auto e: m.edges()) {
        gmds::Edge edg = m.get<gmds::Edge>(e);
        std::vector<TCellID> v_face;
        edg.delegateGetFaceIDs(v_face);
        bool run = true;
        for (auto it: v_face) {
            if (face_boundary->value(it) == 1) {
                edge_boundary->set(e, 1);
                break;
            }
        }
    }
    // Allocate and initializing values
    // Assign boolean value to face (In solution or not) 0 no 1 yes
    this->face_in_solution = m.newVariable<int, gmds::GMDS_FACE>("face_in_solution");
    face_in_solution->setValuesTo(0);
    // Assign face_id for constraint selection
    this->face_id = m.newVariable<int, gmds::GMDS_FACE>("face_id");
    face_id->setValuesTo(0);
    // Assign face_candidates
    this->face_candidates = m.newVariable<int, gmds::GMDS_FACE>("face_candidates");
    face_candidates->setValuesTo(0);
    // Assign boolean value to face (Constrained Strong or not) 0 no 1 yes
    this->face_constraint = m.newVariable<int, gmds::GMDS_FACE>("face_constraint");
    face_constraint->setValuesTo(0);
    // Set Face ID
    std::cout << "Faces" << std::endl;
    for (auto i: m.faces()) {
        this->face_id->set(i, i);
    }

    this->front = m.newVariable<int, gmds::GMDS_FACE>("front");
    this->front->setValuesTo(0);
    std::cout << "#Front " << front->getNbValues() << std::endl;
    std::cout << "#Faces " << m.getNbFaces() << std::endl;
    int count = 0;
    // Ensure that our indexing is fine
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        assert(*it == count);
        count++;
    }
    this->pheromone = m.newVariable<double, gmds::GMDS_FACE>("pheromone");
    this->pheromone->setValuesTo(0.0f);
    this->heuristic = m.newVariable<double, gmds::GMDS_FACE>("heuristic");
    this->heuristic->setValuesTo(1.0f);
    this->probabilities = m.newVariable<double, gmds::GMDS_FACE>("probabilities");
    probabilities->setValuesTo(0.0f);
    /* initialize random seed: */
    this->classification_edge = m.newVariable<int, gmds::GMDS_EDGE>("classification_edge");
    this->classification_face = m.newVariable<int, gmds::GMDS_FACE>("classification_face");
    this->classification_edge->setValuesTo(0);
    this->classification_face->setValuesTo(0);
    this->bnd_solution = m.newVariable<int,gmds::GMDS_EDGE>("bnd_solution");
    this->bnd_solution->setValuesTo(0);

    do_classification_edge();
    do_classification_face();
    srand(time(NULL));
}

std::vector<gmds::TCellID> Env::getEdgesOfFront(std::vector<gmds::TCellID> front, std::vector<char> solution) {
    // Ajouter selection
    std::set<TCellID> edges_set;
    int count;
    for (auto face: front) {
        gmds::Face f = this->m.get<gmds::Face>(face);
        std::vector<gmds::Edge> v_edges;
        std::vector<TCellID> incident_faces;
        f.delegateGet(v_edges);
        //f.getIDs<gmds::Edge>(v_edges);
        for (auto e: v_edges) {
            count = 0;
            e.getIDs<gmds::Face>(incident_faces);
            for (auto face: incident_faces) {
                //if (face_candidates->value(face) != 0) {
                //if (this->front->value(face) > FaceStatus::Candidate) {
                if (solution[face] == 1) {
                    count++;
                }
            }
            if (count > 1) { //! Edge the 'border' of the selecion only
                break;
            }
            if (!this->face_boundary->value(face) && incident_faces.size() == 3) {        //! Case(1)
                break;
            } else if (this->face_boundary->value(face) && incident_faces.size() == 2) {    //! Case(2)
                break;
            }
            edges_set.insert(e.id());
        }
    }
    std::vector<TCellID> edges(edges_set.begin(), edges_set.end());
    return edges;
}

std::vector<gmds::TCellID> Env::getFacesCandidates(std::vector<gmds::TCellID> edges, std::vector<char> solution) {
    std::vector<char> cand(m.getNbFaces(), 0);
    std::set<TCellID> candidates_set;
    int valid = 0;
    for (auto e: edges) { // Edges of front
        std::vector<TCellID> v_face;
        gmds::Edge edg = m.get<gmds::Edge>(e);
        edg.getIDs<gmds::Face>(v_face);
        for (auto f: v_face) { // Candidates
            /*if (solution[f] == 0){
                candidates_set.insert(f);
                cand[f] = 1;
            }*/
            // Remove that thing pls
            gmds::Face face = m.get<gmds::Face>(f);
            std::vector<TCellID> edge_bis;
            face.getIDs<gmds::Edge>(edge_bis);
            for (auto ebis: edge_bis) { // Edges of candidates
                gmds::Edge edge_bis2 = m.get<gmds::Edge>(ebis);

                std::vector<TCellID> v_face_bis;
                edge_bis2.getIDs<gmds::Face>(v_face_bis);
                int counter = 0;
                for (auto fbis: v_face_bis) { // Faces of candidate edges
                    if (solution[fbis] == 1) {
                        counter++;
                    }
                }
                if (counter >= 2) {
                    break;
                }
                if (solution[f] == 0) {
                    candidates_set.insert(f);
                }
            }
        }
    }
    std::vector<TCellID> candidates(candidates_set.begin(), candidates_set.end());
    return candidates;
}

void Env::updateFront(std::vector<TCellID> Front) {
    for (auto f: Front) {
        this->front->set(f, FaceStatus::Front);
    }
}

void Env::updateFrontInitial(std::vector<TCellID> FrontInitial) {
    for (auto f: FrontInitial) {
        this->front->set(f, FaceStatus::Initial);
    }
}

void Env::updateCandidates(std::vector<TCellID> candidates) {
    this->front->setValuesTo(0);
    for (auto c: candidates) {
        this->front->set(c, FaceStatus::Candidate);
    }
}

std::list<gmds::TCellID> Env::getFirstEdgeFront(std::vector<gmds::TCellID> front_initial, std::vector<char> solution) {
    std::list<gmds::TCellID> edge_front;
    for (auto f: front_initial) {
        gmds::Face face = m.get<gmds::Face>(f);
        std::vector<TCellID> edges;
        face.getIDs<gmds::Edge>(edges);
        for (auto e: edges) {
            gmds::Edge cur_edge = m.get<gmds::Edge>(e);
            std::vector<TCellID> faces;
        }
    }

}

std::vector<char> Env::runAnts(std::vector<gmds::TCellID> front_initial) {
    std::vector<char> solution(m.getNbFaces(), false);
    for (auto x: front_initial) {
        solution[x] = 1;
    }
    this->face_constraint_ids = front_initial;
    std::list<gmds::TCellID> edge_front = getFirstEdgeFront(front_initial, solution);//(src.begin(), src.end());
}

void Env::setup_ant_colony_params(int alpha, int beta, float p, float tmin, float tmax, int nb_ants, int iter_global,int n_run) {
    this->alpha = alpha;
    this->beta = beta;
    this->p = p;
    this->tmin = tmin;
    this->tmax = tmax;
    this->nb_ants = nb_ants;
    this->iter_global = iter_global;
    this->n_run = n_run;
}

void Env::resetVariables() {
    face_candidates->setValuesTo(0);
    face_constraint->setValuesTo(0);
    //face_in_solution->setValuesTo(0);
    pheromone->setValuesTo(0.0f);
    heuristic->setValuesTo(0.0f);
    /*
    front->setValuesTo(0);
    for (auto x: face_constraint_ids) {
        front->set(x, FaceStatus::Initial);
    }
    */
}

void Env::show_tcell_vector(std::vector<gmds::TCellID> v) {
    if (this->verbosity == Verbosity::Yes) {
        for (auto e: v) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<TCellID> Env::getEdgesOpposite(TCellID face) {
    std::vector<TCellID> edges_opposite;
    std::vector<TCellID> v_nodes;
    gmds::Face fa = m.get<gmds::Face>(face);
    fa.getIDs<gmds::Node>(v_nodes);
    v_nodes.push_back(v_nodes[0]);
    for (auto it = v_nodes.begin(); it != v_nodes.end() - 1; ++it) {
        std::cout << *it << std::endl;
        std::vector<TCellID> v_edges;
        gmds::Node no = m.get<gmds::Node>(*it);
        no.getIDs<gmds::Edge>(v_edges);
        for (auto x: v_edges) {
            std::vector<TCellID> v_edge_node;
            gmds::Edge e = m.get<gmds::Edge>(x);
            e.getIDs<gmds::Node>(v_edge_node);
            std::cout << "edges  " << v_edge_node.size() << std::endl;
            if (v_edge_node[0] == *it && v_edge_node[1] == *std::next(it) ||
                v_edge_node[1] == *it && v_edge_node[0] == *std::next(it)) {
                edges_opposite.push_back(x);
                break;
            }
        }
    }
    std::cout << "Unordered edges " << std::endl;
    std::vector<TCellID> v_edges_bis;
    gmds::Face facex = m.get<gmds::Face>(face);

    facex.getIDs<gmds::Edge>(v_edges_bis);
    for (auto e: v_edges_bis) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
    std::cout << "Edges edges_opposite " << std::endl;
    for (auto e: edges_opposite) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
    return edges_opposite;
}


void Env::show_normals(std::vector<TCellID> faces) {
    std::cout << "Normals" << std::endl;
    for (auto f: faces) {
        gmds::Face fac = m.get<gmds::Face>(f);
        math::Vector3d normal = fac.normal();
        std::cout << fac.normal() << std::endl;
    }
}

void Env::execute() {
    std::ofstream results_file, all_solutions;
    results_file.open("results_file.txt");
    results_file << "#Parameters " << "p " << this->p << " tmin " << this->tmin << " tmax " << this->tmax << " iteration "
           << this->iter_global << " nb_ants " << this->nb_ants << " n_run " << this->n_run << std::endl;
    results_file << "Run;Iteration;Solution;Quality;Validity;NbEdgeTurn;NbVertexTurn;FaceSelected" << std::endl;

    all_solutions.open("all_solutions.txt");
    all_solutions << "Run;" << "Iter;" << "AntNb;" << "Solution;" << std::endl;

    std::vector<std::vector<char>> best_solution_history;
    std::vector<float> best_solution_quality_history;
    std::vector<int> solution_validity_history;
    /*std::set<gmds::TCellID> ini_faces;
    for(auto s: initial_faces ){
        ini_faces.insert(s.first);
        ini_faces.insert(s.second);
    }
    std::vector<TCellID> initial_faces_vec(ini_faces.begin(), ini_faces.end());*/
    std::ofstream myfile;
    myfile.open("solution_logs.txt");
    std::vector<std::vector<char>> solution;
    std::vector<int> pheromone_add(getNbFaces(), 0);
    std::vector<char> s_best(getNbFaces(), 0);
    float quality_s_best;
    std::vector<int> ids_solutions;
    std::vector<float> pheromones(getNbFaces(), tmax);
    std::vector<std::size_t> valid_solution;
    std::vector<float> qualities;
    std::vector<char> sk_best_solution;
    float sk_best_quality;
    int counter = 0;
    quality_s_best = 0.0f, sk_best_quality = 0.0f;
    std::size_t nb_solution = 0;
    std::vector<char> current_solution;
    for(std::size_t r=0 ; r < this->n_run ; r++ ) {
        for (int i = 0; i < iter_global; i++) {
            // Generate Solutions
            solution.clear();
            valid_solution.clear();
            std::cout << "#Iter " << i << ' ' << this->iter_global << std::endl;
            myfile << "#Iter " << i << std::endl;
            for (int j = 0; j < nb_ants; j++) {
                //solution.push_back(this->run_ant(pheromones)); // Build solution
                current_solution = this->run_ant_bis(pheromones);
                solution.push_back(current_solution); // Build solution
                all_solutions << r << ";" <<  i << ";" << j << ";" << solution_to_string(current_solution) << ";" << std::endl; ;
                ids_solutions.push_back(counter);
                std::cout << "#Solution " << counter << std::endl;
                counter++;
            }
            if (this->verbosity == Verbosity::Yes) {
                for (auto s: solution) {
                    for (auto b: s) {
                        std::cout << b;
                    }
                    std::cout << std::endl;
                }
            }
            int nb_valid_solution = 0;
            for (auto it = solution.begin(); it != solution.end(); ++it) {
                std::size_t size = it->size();
                //if (is_solution_valid(*it)) {
                std::size_t current = std::distance(solution.begin(), it);
                valid_solution.push_back(current);
                if (is_solution_valid(*it)) {
                    nb_valid_solution++;
                }
                //myfile << solution_to_string(solution[current]) << std::endl ;
                //}
            }
            //myfile << "#\tvalid solution " << valid_solution.size() << " / " << solution.size() << std::endl;
            myfile << "#\tvalid solution " << nb_valid_solution << " / " << solution.size() << std::endl;
            if (!valid_solution.empty()) {
                int elite_k = 10;
                // Recuperation des k meilleurs uniques solutions dans un vecteur unique_solution
                std::vector<float> qualities(valid_solution.size(), 0.0f);
                for (auto i_val = 0; i_val < valid_solution.size(); i_val++) {
                    qualities[i_val] = solution_quality(solution[i_val]);
                }
                std::cout << "Valid " << valid_solution.size() << " Total sol " << solution.size() << std::endl;
                assert(valid_solution.size() == solution.size());
                auto x = get_k_best_solutions(10, solution, qualities);
                float quality = 0.0f;
                sk_best_quality = 0.0f;
                std::size_t sk_best_index;
                for (auto it = x.begin(); it != x.end(); ++it) {
                    quality = it->second;
                    if (quality >= sk_best_quality) {
                        sk_best_quality = quality;
                        sk_best_solution = it->first;
                        sk_best_index = std::distance(x.begin(), it);
                    }

                    myfile << "\t\t" << solution_to_string(it->first) << " " << quality << " "
                           << best_quality(it->first) << std::endl;
                }
                // Update pheromon with the new best if exist
                if ((sk_best_quality >= quality_s_best) &&
                    !std::equal(sk_best_solution.begin(), sk_best_solution.end(), s_best.begin())) {
                    quality_s_best = sk_best_quality;
                    s_best = sk_best_solution;
                    // Remove the new best from the set
                    std::swap(x[sk_best_index], x[x.size() - 1]);
                    x.pop_back();
                }
                // Update pheromons with others
                for (auto s: x) {
                    for (auto it = sk_best_solution.begin(); it != sk_best_solution.end(); ++it) {
                        if (*it == 1) {
                            pheromone_add[std::distance(sk_best_solution.begin(), it)] =
                                    1.0 / (1.0 + quality_s_best - sk_best_quality);
                        }
                    }
                }
            }
            best_solution_history.push_back(s_best);
            best_solution_quality_history.push_back(quality_s_best);
            solution_validity_history.push_back(valid_solution.size());
            //! Deposit
            std::string buff;
            for (auto x: pheromones) {
                buff += std::to_string(x) + " ";
            }
            buff += "\n";
            for (std::size_t i = 0; i < pheromones.size(); i++) {
                pheromones[i] = pheromones[i] * p + pheromone_add[i];
                if (pheromones[i] > tmax) {
                    pheromones[i] = tmax;
                } else if (pheromones[i] < tmin) {
                    pheromones[i] = tmin;
                }
                buff += std::to_string(pheromones[i]) + " ";
            }
            myfile << buff + "\n";
            updatePheromones(pheromones);
            //writeSolution(s_best, "best" + format_number(i, 5) + ".vtk", pheromones, {}, {});
            std::fill(pheromone_add.begin(), pheromone_add.end(), 0.0f);
        }

        myfile << "#Best solution\t" << solution_to_string(s_best) << " Q= " << solution_quality(s_best)
               << " egde_turn "
               << nb_edge_turn(s_best) << std::endl;
        myfile.close();
        set_write(true);
        writeSolution(s_best, std::to_string(r)+"best_solution.vtk", {}, {}, {});
        set_write(false);
        std::cout << "#Best solution\t" << solution_to_string(s_best) << " Q= " << solution_quality(s_best) << " valid "
                  << is_solution_valid(s_best) << std::endl;
        for (std::size_t i = 0; i < best_solution_history.size(); i++) {
            results_file << r << ";" << i << ";" << solution_to_string(best_solution_history[i]) << ";"
                   << best_solution_quality_history[i]
                   << ";" << solution_validity_history[i] << ";" << nb_edge_turn(best_solution_history[i]) << ";"
                   << nb_vertex_turn(best_solution_history[i]) << ";"
                   << std::count(best_solution_history[i].begin(), best_solution_history[i].end(), 1) << std::endl;
        }
        best_solution_history.clear();
        solution_validity_history.clear();
    }
    myfile.close();
    results_file.close();
    all_solutions.close();
}


std::vector<char> Env::run_ant(std::vector<float> pheromones) {
    std::ofstream myfile;
    myfile.open("logs.txt");
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    char *dt = ctime(&now);
    std::vector<std::pair<TCellID, TCellID>> edge_front(this->initial_front);
    myfile << dt << "\n";

    std::vector<char> solution = constructInitialSolution(this->initial_faces);
    std::vector<char> processed_edge(m.getNbEdges());
    std::vector<TCellID> init_faces = solution_to_faces(solution);

    std::fill(processed_edge.begin(), processed_edge.end(), 0);
    std::cout << "run_ant" << std::endl;
    std::vector<TCellID> history_face;
    std::vector<TCellID> history_edge;
    init_disable_faces();

    int count = 0;
    bool run = true;
    while (!edge_front.empty() && run ) {
        myfile << "Step " << count << " edges_in_front: ";
        for (auto e: edge_front) {
            myfile << e.first << ' ';
        }
        myfile << '\n';
        std::vector<std::pair<TCellID, TCellID>> candidates;
        std::map<TCellID ,std::vector<TCellID>> candidates_map;
        candidates_map = getCandidate_v3(edge_front,solution);
        myfile << "#candidates " << candidates_map.size() << std::endl;
        /*for (auto c: candidates) {
            myfile << c.first << ' ';
        }*/
        for (auto c: candidates_map) {
            myfile << c.first << "\n\t";
            for(auto e : c.second) {
                myfile << e << ' ';
            }
            myfile << "\n" ;
        }
        myfile << std::endl;
        myfile << "S: ";
        for (auto s: solution) {
            myfile << std::to_string(s);
        }
        myfile << "\n";
        std::cout << "here" << std::endl;
        std::vector<std::pair<TCellID, TCellID>> new_edges;
        std::vector<std::pair<TCellID, float>> proba;
        if (!candidates_map.empty()) {
            //std::pair<TCellID,TCellID> selected = candidates[0],tmp;
            std::pair<TCellID, TCellID> selected, tmp;
            if (count == 0) {
                myfile << "Random selection " ;
                selected = random_selection(candidates_map);
            } else {

                for(auto f_es : candidates_map){
                    candidates.push_back(std::make_pair(f_es.first,f_es.second[0]));
                }
                selected = selectCandidate(candidates, solution, pheromones, proba);
                myfile << "Face:Proba ," ;
                for(auto f_p : proba){
                    myfile << f_p.first << ":" << f_p.second << ", " ;
                }
                myfile << "\n" ;
            }
            //selectCandidate_bis(candidates,solution);
            myfile << "selected " << selected.first << " from edge " << selected.second << "\n";
            solution[selected.first] = 1;
            history_face.push_back(selected.first);
            history_edge.push_back(selected.second);
            auto it = std::find_if(edge_front.begin(), edge_front.end(),
                                   [&selected](const std::pair<TCellID, TCellID> &element) {
                                       return element.first == selected.second;
                                   });


            myfile << "removed edge" << " " ;
            disable_faces_of_selected(candidates,selected);
            // Remove the edge of the candidate processed
            auto edge_to_remove = candidates_map[selected.first];
            for( auto rmv : edge_to_remove ){
                auto it = std::find_if( edge_front.begin(), edge_front.end(),
                                        [&rmv](const std::pair<TCellID , TCellID>& element){ return element.first == rmv;} );
                if( it != edge_front.end()) {
                    tmp = edge_front.back();
                    std::swap(edge_front[std::distance(edge_front.begin(), it)],edge_front[edge_front.size()-1]);
                    myfile << it->first << " " ;
                    edge_front.pop_back();
                    processed_edge[it->first] = 1;
                }else{
                    assert(1==2);
                }
            }
            //new_edges = getEdges_v2(selected.first, processed_edge, solution);
            new_edges = getEdges_v3(selected.first,selected.second,processed_edge);
            myfile << "\nedges to add: ";
            for (auto e: new_edges) {
                myfile << e.first << ' ';
            }
            myfile << '\n';
            edge_front.insert(edge_front.end(), new_edges.begin(), new_edges.end());

        } else {
            run = false;
        }
        count++;
        std::cout << proba.size() << std::endl;

        /*
        update_probabilities(proba);
        updatePheromones(pheromones);
        updateCandidates(cdt);
        updateFront(faces_history);
        writeSolution(solution, "exec" + format_number(count, 5) + ".vtk", {});
        std::cout << "#S valid, run, edge_front " << !is_solution_valid(solution) << ", " << run << ", "
                  << edge_front.size() << std::endl;
        myfile << "#S valid, run, edge_front " << !is_solution_valid(solution) << ", " << run << ", "
               << edge_front.size() << std::endl;
        */
        std::cout << "loop done" << std::endl;
    }

    // Last write

    //updateCandidates(cdt);
    updateFront(history_face);
    assert(this->initial_front.size()!=0);
    writeSolution(solution, "exec" + format_number(count, 5) + ".vtk", {},history_face,history_edge);

    myfile << " Nb_edge_turn " << nb_edge_turn(solution) << " nb edg " << m.getNbEdges() << std::endl;
    myfile << " Nb_vtx_turn " << nb_vertex_turn(solution) << std::endl;
    myfile << " Selection validity " << is_selection_valid(solution) << std::endl;
    myfile << " Solution validity " << is_solution_valid(solution) << std::endl;
    myfile << " Quality " << solution_quality(solution) << "/" << best_quality(solution) << std::endl;
    for (auto s: solution) {
        myfile << std::to_string(s);
    }
    myfile.close();
    return solution;
}

std::pair<gmds::TCellID, gmds::TCellID>
Env::selectCandidate(std::vector<std::pair<gmds::TCellID, gmds::TCellID>> candidates, std::vector<char> solution,
                         std::vector<float> pheromones, std::vector<std::pair<TCellID, float>> &face_proba) {
    assert(candidates.size() != 0);
    std::cout << "#Candidates " << candidates.size() << " ";
    for (auto c: candidates) {
        std::cout << c.first << " ";
    }
    float sum_all = 0.0f;
    //! Compute total pheromones
    float sum_pheromones = 0.0f;
    for (auto c: candidates) {
        sum_pheromones += pheromones[c.first];
    }
    std::cout << "sum pheromones " << sum_pheromones << " for " << candidates.size() << " candidates" << std::endl;
    std::accumulate(pheromones.begin(), pheromones.end(), 0);
    std::vector<std::pair<float, float>> intervals;
    intervals.reserve(candidates.size());
    std::vector<float> proba;
    std::vector<float> heuristic;
    proba.reserve(candidates.size());
    float prev = 0.0f, cur;

    // Pheromone factor
    float pheromones_sum = 0.0f;
    for (auto c: candidates) {
        pheromones_sum += pheromones[c.first];
    }
    // heuristicfactor
    float heuristic_sum = 0.0f;
    for (auto c: candidates) {
        if (same_orientation(c.first, c.second, solution)) {
            heuristic.push_back(2.0);
        } else {
            heuristic.push_back(1.0);
        }
    }
    for (std::size_t i = 0; i < candidates.size(); i++) {
        sum_all += pow(pheromones[candidates[i].first], alpha) * pow(heuristic[i], beta);
    }
    for (auto h: heuristic) {
        heuristic_sum += pow(h, beta);
    }

    // Compute probability
    float previous = 0.0f, current;
    for (std::size_t i = 0; i < candidates.size(); i++) {
        float ph = pheromones[candidates[i].first];
        float h = heuristic[i];
        current = (pow(ph, alpha) * pow(h, beta)) / sum_all;
        assert(0.0 <= current && current <= 1.0);
        proba.push_back(current);
        face_proba.push_back(std::make_pair(candidates[i].first, current));
        previous = current;
    }
    assert(heuristic.size() == candidates.size());
    std::cout << "PB1";
    float sum = 0.f;
    previous = 0.f, current = 0.0f;
    for (auto it: proba) {
        current = it;
        sum += current;
        intervals.push_back(std::make_pair(previous, sum));
        previous = sum;
    }
    std::cout << "Total " << sum << std::endl;
    std::cout << "PB2";
    for (auto it: intervals) {
        std::cout << '(' << it.first << ", " << it.second << ") ";
    }
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); // Float in  [0,1]
    std::cout << "Interval " << r << std::endl;
    assert(0.f <= r && r <= 1.0f);
    std::pair<gmds::TCellID, gmds::TCellID> selected;
    bool found = false;
    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
        if (((*it).first <= r) && (r <= (*it).second)) {
            selected = candidates[std::distance(intervals.begin(), it)];
            std::cout << "selected " << selected.first << std::endl;
            found = true;
            break;
        }
    }
    assert(found == true);
    return selected;
}