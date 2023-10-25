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
    if (this->write_vtk) {
        IGMeshIOService ios(&(this->m));
        VTKWriter writer(&ios);
        writer.setCellOptions(gmds::N | gmds::E | gmds::F | gmds::R);
        writer.setDataOptions(gmds::N | gmds::E | gmds::F | gmds::R);
        writer.write(
                fname); // paraview ctrl + space (Shrink) pour avoir uniquement les hex (enlever les F dans l'écriture)
    }
}

void Env::writeSolution(std::vector<char> solution, std::string fname,
                        std::vector<float> pheromones, std::vector<TCellID> history_faces,
                        std::vector<TCellID> history_edges) {


    Variable<int> *history_faces_visu = m.newVariable<int, gmds::GMDS_FACE>("history_faces");
    Variable<int> *history_edges_visu = m.newVariable<int, gmds::GMDS_EDGE>("history_edges");
    history_faces_visu->setValuesTo(0);
    history_edges_visu->setValuesTo(0);
    int cpt = 0;
    for (auto x: history_faces) {
        history_faces_visu->set(x, cpt);
        cpt++;
    }
    cpt = 0;
    for (auto x: history_edges) {
        history_edges_visu->set(x, cpt);
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
    m.deleteVariable(gmds::GMDS_FACE, "history_faces");
    m.deleteVariable(gmds::GMDS_EDGE, "history_edges");
}

Env::Env(std::string fname) : m(gmds::MeshModel(
        DIM3 | F | R | R2N | R2F | F2R | E | N | F2N | E2N |
        F2E | E2F | N2F | N2E)) {
    std::string vtk_file = fname;
    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::E | gmds::F | R);
    vtkReader.setDataOptions(gmds::N | gmds::E | gmds::F | R);
    vtkReader.read(vtk_file);

    std::cout << "#Nodes " << m.getNbFaces() << std::endl;
    std::cout << "Filename " << fname << std::endl;
    std::cout << "Nodes " << m.getNbNodes() << std::endl;
    std::cout << "Regions " << m.getNbRegions() << std::endl;
    std::cout << "Edges " << m.getNbEdges() << std::endl;
    std::cout << "Faces " << m.getNbFaces() << std::endl;
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        Face f = m.get<Face>(*it);
        std::cout << "(" << f.get<Region>().size() << ", " << f.get<Edge>().size() << ") ";
    }
    std::cout << std::endl;

    gmds::MeshDoctor doc(&m);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();
    std::size_t mcpt = 0;
    std::cout << "#Nodes " << m.getNbFaces() << std::endl;
    std::cout << "Filename " << fname << std::endl;
    std::cout << "Nodes " << m.getNbNodes() << std::endl;
    std::cout << "Regions " << m.getNbRegions() << std::endl;
    std::cout << "Edges " << m.getNbEdges() << std::endl;
    std::cout << "Faces " << m.getNbFaces() << std::endl;
    //std::cout << "Hexaedra " << m.getNbHexahedra() << std::endl;
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        Face f = m.get<Face>(*it);
        std::cout << "(" << f.get<Region>().size() << ", " << f.get<Edge>().size() << ") ";
    }
    std::cout << std::endl;

    // Compute
    if (m.hasVariable(gmds::GMDS_EDGE, "edge_boundary")) {
        this->edge_boundary = m.getVariable<int, gmds::GMDS_EDGE>("edge_boundary");
    } else {
        this->edge_boundary = m.newVariable<int, gmds::GMDS_EDGE>("edge_boundary");
    }
    if (m.hasVariable(gmds::GMDS_FACE, "face_boundary")) {
        this->face_boundary = m.getVariable<int, gmds::GMDS_FACE>("face_boundary");
    } else {
        this->face_boundary = m.newVariable<int, gmds::GMDS_FACE>("face_boundary");
    }
    // Get
    this->face_in_solution = m.getVariable<int, gmds::GMDS_FACE>("face_in_solution");
    // useless
    if (m.hasVariable(gmds::GMDS_FACE, "face_id")) {
        this->face_id = m.getVariable<int, gmds::GMDS_FACE>("face_id");
    } else {
        this->face_id = m.newVariable<int, gmds::GMDS_FACE>("face_id");
    }

    //auto blblbl = m.getVariable<int, gmds::GMDS_FACE>("fblblblbl");
    //! ^--- make the program crash if variable not found whereas it is not rly a problem if we can compute it

    //this->face_constraint = m.getVariable<int, gmds::GMDS_FACE>("face_constraint");
    //this->front = m.getVariable<int, gmds::GMDS_FACE>("front");

    if (m.hasVariable(gmds::GMDS_EDGE, "hex_number_edge")) {
        this->hex_number_edge = m.getVariable<int, gmds::GMDS_EDGE>("hex_number_edge");
    } else {
        assert(1 == 0);
    }
    if (m.hasVariable(gmds::GMDS_FACE, "hex_number_face")) {
        this->hex_number_face = m.getVariable<int, gmds::GMDS_FACE>("hex_number_face");
    } else {
        assert(1 == 0);
    }
    if (m.hasVariable(gmds::GMDS_EDGE, "classification_edge")) {
        this->classification_edge = m.getVariable<int, gmds::GMDS_EDGE>("classification_edge");
    } else {
        this->classification_edge = m.newVariable<int, gmds::GMDS_EDGE>("classification_edge");
        this->classification_edge->setValuesTo(0);
        do_classification_edge();
    }
    if (m.hasVariable(gmds::GMDS_FACE, "classification_face")) {
        this->classification_face = m.getVariable<int, gmds::GMDS_FACE>("classification_face");
    } else {
        this->classification_face = m.newVariable<int, gmds::GMDS_FACE>("classification_face");
        this->classification_face->setValuesTo(0);
        do_classification_face();
    }
    this->bnd_solution = m.newVariable<int, gmds::GMDS_EDGE>("bnd_solution");
    this->bnd_solution->setValuesTo(0);

    if (m.hasVariable(gmds::GMDS_FACE, "pheromone")) {
        this->pheromone = m.getVariable<double, gmds::GMDS_FACE>("pheromone");
    } else {
        this->pheromone = m.newVariable<double, gmds::GMDS_FACE>("pheromone");
        pheromone->setValuesTo(0.0f);
    }

    if (m.hasVariable(gmds::GMDS_FACE, "heuristic")) {
        this->heuristic = m.getVariable<double, gmds::GMDS_FACE>("heuristic");
    } else {
        this->heuristic = m.newVariable<double, gmds::GMDS_FACE>("heuristic");
        heuristic->setValuesTo(0.0f);
    }

    if (m.hasVariable(gmds::GMDS_FACE, "face_candidates")) {
        this->face_candidates = m.getVariable<int, gmds::GMDS_FACE>("face_candidates");
    } else {
        this->face_candidates = m.newVariable<int, gmds::GMDS_FACE>("face_candidates");
        this->face_candidates->setValuesTo(0.0f);
    }
    if (m.hasVariable(gmds::GMDS_FACE, "face_constraint")) {
        this->face_constraint = m.getVariable<int, gmds::GMDS_FACE>("face_constraint");
    } else {
        this->face_constraint = m.newVariable<int, gmds::GMDS_FACE>("face_constraint");
        this->face_constraint->setValuesTo(0.0f);
    }
    // Check stuff
    // Number hex for
    // faces
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        assert(this->hex_number_face->value(*it) > 0);
    }
    //  edges
    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        assert(this->hex_number_edge->value(*it) > 0);
    }
    // Check connectivity
    // E2N
    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        Edge e = m.get<Edge>(*it);
        std::cout << *it << std::endl;
        assert(e.nbNodes() > 0);
    }
    // F2E
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        Face f = m.get<Face>(*it);
        std::cout << *it << std::endl;
        assert(f.nbEdges() > 0);
    }
    // E2F
    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        Edge e = m.get<Edge>(*it);
        std::cout << *it << " " << e.getIDs<Face>().size() << std::endl;
        assert(e.nbFaces() > 0);
    }

    for (auto i = 0; i < this->face_in_solution->getNbValues(); i++) {
        if (this->face_in_solution->value(i) != 0) {
            this->face_in_solution->set(i, 3);
        }
        std::cout << this->face_in_solution->value(i) << " ";
    }

    writeVTK("test.vtk");
    srand(time(NULL));
}

void Env::do_classification_edge(){
    int s_region = m.getNbRegions();
    for (auto e: m.edges()) {
        Edge edge = m.get<Edge>(e);
        std::vector<TCellID> regions;
        //edge.getIDs<Region>(regions);
        //assert(edge.nbRegions() == regions.size());
        //int size_e = m.edges().getNbElements();
        //int size_e_r = edge.nbRegions();
        //int size_e_r0 = regions.size();
        int classification;
        switch (this->hex_number_edge->value(e)) {
            case 3:
                classification = (int) edge_classification::concave;    // 0
                break;
            case 1:
                classification = (int) edge_classification::convex;     // 1
                break;
            case 2:
                classification = (int) edge_classification::surface;    // 2
                break;
            case 4:
                classification = (int) edge_classification::volume;     // 3
                break;
            default:
                assert(1 == 0);
        }
        classification_edge->set(e, classification);
    }
}

void Env::do_classification_face(){
    for (auto f: m.faces()) {
        Face face = m.get<Face>(f);
        //std::vector<TCellID> regions;
        //face.getIDs<Region>(regions);
        int classification;
        //std::cout << "here " << regions.size() << " " << this->hex_number_face->value(f) << " " << face.nbRegions() << std::endl;
        //assert(regions.size() == this->hex_number_face->value(f));
        switch (this->hex_number_face->value(f)){
            case 1:
                classification = (int) face_classification::surface;
                break;
            case 2:
                classification = (int) face_classification::volume;
                break;
            default:
                assert(1 == 0);
        }
        //std::cout << f << " " << face.nbRegions() << " " << regions.size() <<  std::endl;
        classification_face->set(f, classification);
    }
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
    this->bnd_solution = m.newVariable<int, gmds::GMDS_EDGE>("bnd_solution");
    this->bnd_solution->setValuesTo(0);

    this->hex_number_edge = m.newVariable<int, gmds::GMDS_EDGE>("hex_number_edge");
    this->hex_number_edge->setValuesTo(0);
    this->hex_number_face = m.newVariable<int, gmds::GMDS_FACE>("hex_number_face");
    this->hex_number_face->setValuesTo(0);

    initialize_hex_number_edge();
    initialize_hex_number_face();

    do_classification_edge();
    do_classification_face();

    srand(time(NULL));

    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        Edge e = m.get<Edge>(*it);
        assert(this->hex_number_edge != 0);
    }
    for (auto it = m.faces().begin(); it != m.faces().end(); ++it) {
        Face f = m.get<Face>(*it);
        assert(this->hex_number_face != 0);
    }
    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        Edge e = m.get<Edge>(*it);
        std::cout << *it << " " << e.getIDs<Face>().size() << std::endl;
        assert(e.nbFaces() > 0);
    }
    for (auto it = m.edges().begin(); it != m.edges().end(); ++it) {
        Edge e = m.get<Edge>(*it);
        std::cout << *it << std::endl;
        assert(e.nbNodes() > 0);
    }
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


void Env::setup_ant_colony_params(int alpha, int beta, float p, float tmin, float tmax, int nb_ants, int iter_global,
                                  int n_run, int nb_elite, float tau_initial) {
    this->alpha = alpha;
    this->beta = beta;
    this->p = p;
    this->tmin = tmin;
    this->tmax = tmax;
    this->nb_ants = nb_ants;
    this->iter_global = iter_global;
    this->n_run = n_run;
    this->nb_elite = nb_elite;
    this->tau_initial = tau_initial;
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
    std::ofstream results_file, all_solutions,myfile,file_csv,pheromones_csv,recap_file;
    recap_file.open("recap_file.txt");
    recap_file << "\np" << this->p << "\ntmin " << this->tmin << "\ntmax " << this->tmax << "\niteration " << this->iter_global << "\nnb_ants " << this->nb_ants << "\nn_run " << this->n_run << "\ntini " << this->tau_initial << "\nn_elite " << this->nb_elite << std::endl;
    recap_file << "\nModel:\n#Hexa " << m.getNbRegions() << "\n#F " << m.getNbFaces() << "\n#E " << m.getNbEdges() << "\n#N " << m.getNbNodes() << std::endl;
    results_file.open("results_file.txt");

    results_file << "#Parameters " << "p " << this->p << " tmin " << this->tmin << " tmax " << this->tmax
                 << " iteration "
                 << this->iter_global << " nb_ants " << this->nb_ants << " n_run " << this->n_run << std::endl;
    results_file
            << "Run;Iteration;BestSolution;BestQuality;Validity;NbEdgeTurn;NbVertexTurn;FaceSelected;Mean;Q1;Q3;Median;"
            << std::endl;

    all_solutions.open("all_solutions.txt");
    all_solutions << "Run;" << "Iter;" << "AntNb;" << "Solution;" << "NewSolution;" << "Quality;" << std::endl;

    std::ofstream discovered_per_iteration;
    discovered_per_iteration.open("discovered_per_iteration.txt");
    discovered_per_iteration << "Run;" << "Iter;" << "NbSolutionDiscovered;" << std::endl;

    std::vector<std::vector<char>> best_solution_history;
    std::vector<float> best_solution_quality_history;
    std::vector<int> solution_validity_history;

    best_solution_history.reserve(this->iter_global * this->nb_ants);
    best_solution_quality_history.reserve(this->iter_global * this->nb_ants);

    myfile.open("solution_logs.txt");
    file_csv.open("run_results.csv");
    file_csv << "#n_run=" << this->n_run << std::endl;
    file_csv << "i;Q1;Q3;Med;Avg;Best;Valid" << std::endl;
    // Dec

    std::set<std::string> solution_set;


    std::vector<std::vector<char>> solution;
    float quality_s_best;
    std::vector<int> ids_solutions;
    std::vector<float> pheromones(getNbFaces(), this->tau_initial);
    std::vector<std::size_t> valid_solution;
    std::vector<float> qualities;
    std::vector<char> sk_best_solution;
    std::vector<char> s_best(getNbFaces(), 0);
    std::vector<float> pheromone_add(getNbFaces(), 0.0);
    float sk_best_quality;
    int counter = 0;
    std::size_t nb_solution = 0;
    std::vector<char> current_solution;
    std::vector<float> q1, q2, mean, median;
    std::vector<float> all_solution;
    /** Stat **/
    std::vector<std::vector<float>> all_qualities(this->iter_global);

    for (std::size_t r = 0; r < this->n_run; r++) {
        quality_s_best = 0.0f, sk_best_quality = 0.0f;
        std::size_t nb_unique_solution = 0;
        std::fill(pheromones.begin(), pheromones.end(), this->tau_initial);
        pheromones_csv.open(std::to_string(r)+"pheromone.csv" );
        pheromones_csv << "iteration;total;added;part;maximum_ph;maximum_ph_add;" << std::endl;
        for (int i = 0; i < iter_global; i++) {
            // Generate Solutions
            solution.clear();
            valid_solution.clear();
            std::cout << "#Iter " << i << ' ' << this->iter_global << std::endl;
            myfile << "#Iter " << i << std::endl;
            for (int j = 0; j < nb_ants; j++) {
                current_solution = this->run_ant_bis(pheromones);
                solution.push_back(current_solution); // Build solution
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
                // Recuperation des k meilleurs uniques solutions dans un vecteur unique_solution
                std::vector<float> qualities(valid_solution.size(), 0.0f);
                for (auto i_val = 0; i_val < valid_solution.size(); i_val++) {
                    qualities[i_val] = solution_quality(solution[i_val]);
                    std::cout << qualities[i_val] << std::endl;
                    std::string str_sol = solution_to_string(solution[i_val]);
                    auto inserted = solution_set.insert(str_sol);
                    if( inserted.second){
                        nb_unique_solution++;
                    }
                    all_solutions << r << ";" << i << ";" << str_sol << ";" << inserted.second << ";" << qualities[i_val] << ";" << std::endl;
                    //assert(qualities[i_val] <= 3.0 && qualities[i_val] >= 0.0);
                }
                discovered_per_iteration << r << ";" << i << ";" << nb_unique_solution << ";" << std::endl ;
                std::cout << "Valid " << valid_solution.size() << " Total sol " << solution.size() << std::endl;
                assert(valid_solution.size() == solution.size());
                auto x = get_k_best_solutions(nb_elite, solution, qualities);
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
                    myfile << "\t\t" << solution_to_string(it->first) << " " << quality << " " << best_quality(it->first) << std::endl;
                }
                // Update pheromon with the new best if exist
                if ((sk_best_quality >= quality_s_best)){
                    //&& !std::equal(sk_best_solution.begin(), sk_best_solution.end(), s_best.begin())) {
                    quality_s_best = sk_best_quality;
                    s_best = sk_best_solution;
                    // Remove the new best from the set
                    //std::swap(x[sk_best_index], x[x.size() - 1]);
                    //x.pop_back();
                }
                // Update pheromons with others
                for (auto s: x){
                    for( auto face_it = s.first.begin() ; face_it != s.first.end() ; ++face_it ){
                        if( *face_it == 1 ){
                            std::size_t index = std::distance(s.first.begin(), face_it);
                            pheromone_add[index] += 1.0 / (1.0 + quality_s_best - s.second);
                        }
                    }
                }
                std::copy(qualities.begin(), qualities.end(), std::back_inserter(all_qualities[i]));
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
            float ph_tot =  std::accumulate(pheromones.begin(),pheromones.end(),0.0);
            float ph_add =  std::accumulate(pheromone_add.begin(),pheromone_add.end(),0.0);
            float maximum_ph_on_face = *std::max_element(pheromones.begin(),pheromones.end());
            float maximum_ph_add = *std::max_element(pheromone_add.begin(),pheromone_add.end());
            pheromones_csv << i << ";" << ph_tot << ";" << ph_add << ";" << ph_tot / ph_add << ";" << maximum_ph_on_face << ";" << maximum_ph_add << ";" << std::endl;
            std::fill(pheromone_add.begin(), pheromone_add.end(), 0.0f);
        }

        myfile << "#Best solution\t" << solution_to_string(s_best) << " Q= " << solution_quality(s_best)
               << " egde_turn "
               << nb_edge_turn(s_best) << std::endl;
        myfile.close();
        set_write(true);
        writeSolution(s_best, std::to_string(r) + "best_solution.vtk", {}, {}, {});
        set_write(false);
        std::cout << "#Best solution\t" << solution_to_string(s_best) << " Q= " << solution_quality(s_best) << " valid "
                  << is_solution_valid(s_best) << std::endl;
        for (std::size_t i = 0; i < best_solution_history.size(); i++) {
            results_file << r << ";" << i << ";" << solution_to_string(best_solution_history[i]) << ";"
                         << best_solution_quality_history[i]
                         << ";" << solution_validity_history[i] << ";" << nb_edge_turn(best_solution_history[i]) << ";"
                         << nb_vertex_turn(best_solution_history[i]) << ";"
                         << std::count(best_solution_history[i].begin(), best_solution_history[i].end(), 1)
                         << std::endl;
        }
        pheromones_csv.close();
        best_solution_history.clear();
        solution_validity_history.clear();
    }
    int count = 0;
    for (auto x: all_qualities) {
        std::sort(x.begin(), x.end());
    }
    std::size_t med, quart1, quart3;
    float avg;
    for (std::size_t x = 0; x < all_qualities.size(); x++) {
        std::sort(all_qualities[x].begin(), all_qualities[x].end());
        med = static_cast<int>((all_qualities[x].size() + 1) / 2);
        quart1 = static_cast<int>((all_qualities[x].size() + 1) / 4);
        quart3 = static_cast<int>(med + quart1);//"i;Q1;Q3;Med;Avg;Best;Valid"
        avg = std::accumulate(all_qualities[x].begin(), all_qualities[x].end(), 0.0) / all_qualities[x].size();
        file_csv << x << ";" << all_qualities[x][quart1] << ";" << all_qualities[x][quart3] << ";"
                 << all_qualities[x][med] << ";" << avg << ";" << all_qualities[x].back() << ";"
                 << all_qualities[x].size() << std::endl; // i Q1 Q3 Med Avg Best
    }
    myfile.close();
    results_file.close();
    all_solutions.close();
    file_csv.close();
    discovered_per_iteration.close();

}

