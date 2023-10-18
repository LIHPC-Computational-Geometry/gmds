#ifndef ENV_H
#define ENV_H
/******************************************************************************/
#include <unordered_map>
/******************************************************************************/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/ig/Mesh.h>
#include <gmds/math/VectorDyn.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <numeric>
#include <iostream>
#include <stdio.h>
#include <iterator>

#include <fstream>
#include <ctime>
#include <cmath>
/******************************************************************************/
#define _INFINITY 10000

std::string format_number(int number, int max_digit);

enum Verbosity {
    No = 0, Yes = 1
};

struct ordering {
    bool operator()(std::pair<int, float> const &a,
                    std::pair<int, float> const &b) {
        return a.second < b.second;
    }
};

struct ordering_greater {
    bool operator()(std::pair<std::string, float> const &a,
                    std::pair<std::string, float> const &b) {
        return a.second > b.second;
    }
};

enum class edge_classification {
    concave = 0, convex = 1, surface = 2, volume = 3,
};

enum class face_classification {
    surface = 2, volume = 3,
};




namespace gmds {

    class Env {
        enum FaceStatus {
            None = 0, Candidate = 1, Front = 2, Initial = 3
        };

    private:
        // Our mesh
        gmds::Mesh m;
        // Dims
        TInt x_n, y_n, z_n;
        int verbosity;
        // Faces constraints
        std::vector<gmds::TCellID> face_constraint_ids;
        // Assign boolean value to face (In solution or not) 0 no 1 yes
        Variable<int> *face_in_solution;
        // Assign face_id for constraint selection
        Variable<int> *face_id;
        // Assign face_candidates
        Variable<int> *face_candidates;
        // Assign boolean value to face (Constrained Strong or not) 0 no 1 yes
        Variable<int> *face_constraint;
        // Assign boolean value to blocs (exists or not) 0 no 1 yes
        //Variable<int> *bloc_exist;
        Variable<int> *var_delete;

        Variable<int> *edge_boundary;

        Variable<int> *face_boundary; // 1 on boundary
        Variable<int> *front;
        Variable<int> *classification_edge;
        Variable<int> *classification_face;

        Variable<double> *pheromone;
        Variable<double> *probabilities;
        Variable<double> *heuristic;
        Variable<int> *hex_number_edge;
        Variable<int> *hex_number_face;


        Variable<int>* bnd_solution;

        std::vector<TCellID> initial_faces;
        std::vector<std::pair<TCellID, TCellID>> initial_front;

        std::vector<char> disable_faces;

        // Experiment params
        int n_run,nb_elite;
        // Ant colony params
        int alpha, beta, nb_ants, iter_global;
        float p, tmin, tmax, tau_initial;
        // Write vtk verbosity
        int write_vtk=1; //! 0 no  write 1 write only one result 2 write all the build step of the selection (high cost)

    public :

        void set_initial_faces(std::vector<TCellID> initial_faces){
            this->initial_faces = initial_faces;
        }

        void set_write(int write){
            assert( 0 <= write <= 2);
            this->write_vtk = write;
        }
        /** **/
        void reset_solution(){
            for(std::size_t i =0; i < face_in_solution->getNbValues() ; i++){
                if(face_in_solution->value(i) != 3){
                    face_in_solution->value(i) = 0;
                }
            }

        }

        int nb_faces_in_solutions(TCellID edge) {
            gmds::Edge edges = m.get<gmds::Edge>(edge);
            std::vector<TCellID> faces = edges.getIDs<gmds::Face>();
            int cpt = 0;
            for (auto f: faces) {
                if (this->face_in_solution->value(f) != 0) {
                    cpt++;
                    std::cout << f << std::endl;
                }
            }
            return cpt;
        }

        void build_bnd_solution(){
            this->bnd_solution->setValuesTo(0);
            for(auto f : m.faces()){
                Face face = m.get<Face>(f);
                std::vector<TCellID> edges = face.getIDs<gmds::Edge>();
                if(this->face_in_solution->value(f) != 0) {
                    for (auto e: edges) {
                        if (edge_to_add(static_cast<face_classification>(classification_face->value(f)),
                                        static_cast<edge_classification>(classification_edge->value(e)))) {
                            int nb = nb_faces_in_solutions(e);
                            this->bnd_solution->set(e, nb);
                            std::cout << nb << std::endl;
                            // 1 bord
                            // 0 externe
                            // 2 in solution
                            // 3 err
                        }
                    }
                }
            }
        }

        std::set<std::pair<TCellID,TCellID>> getCandidates(){
            std::set<std::pair<TCellID,TCellID>> candidates;
            for(auto e: m.edges()){
                if(this->bnd_solution->value(e) == 1){
                    Edge edge = m.get<Edge>(e);
                    std::vector<TCellID> faces;
                    edge.getIDs<gmds::Face>(faces);
                    for(auto f : faces){
                        if( this->face_in_solution->value(f) == 0 && edge_of_faces_solution(f) ) {
                            candidates.insert(std::make_pair(f,e));
                        }
                    }
                }
            }
            return candidates;
        }

        void initialize_ant_bis(std::vector<TCellID> faces){
            assert(!faces.empty());
            if(!faces.empty()) {
                for (auto f: faces) {
                    this->face_in_solution->set(f, 3);
                }
            }

        }

        void initialize_ant_bis(std::vector<std::pair<TCellID,TCellID>> faces){
            this->face_in_solution->set(faces.front().first,3);
            for(auto f = std::next(faces.begin());  f != faces.end() ; ++f){
                this->face_in_solution->set(f->second,3);
            }
        }

        TCellID select_random(){
            int min = 0, max = m.getNbFaces() - 1;
            int randNum = rand() % (max - min + 1) + min;
            return randNum;
        }


        std::vector<char> run_ant_bis(std::vector<float> pheromones) {
            reset_solution();
            std::vector<char> solution;
            std::set<std::pair<TCellID, TCellID>> candidates;
            TCellID selected;
            int i=0;
            std::vector<std::pair<TCellID, float>> face_proba;
            selected = select_random();
            if( this->face_in_solution->value(selected) != 3 ) {
                this->face_in_solution->set(selected, 1);
            }
            build_bnd_solution();
            candidates = getCandidates();
            for(auto c : candidates){
                this->face_candidates->set(c.first,1);
            }
            writeVTK(format_number(i, 5)+"test_run_ant_bis.vtk");
            this->face_candidates->setValuesTo(0);
            std::cout << "Candidates" << candidates.size() << std::endl;
            while (!candidates.empty()){
                selected = select_Candidate(candidates, pheromones, face_proba);
                this->face_in_solution->set(selected,1);
                build_bnd_solution();
                candidates = getCandidates();
                i++;
                for(auto c : candidates){
                    this->face_candidates->set(c.first,1);
                }
                writeVTK(format_number(i, 5)+"test_run_ant_bis.vtk");
                this->face_candidates->setValuesTo(0);
            }
            for(std::size_t i = 0; i < face_in_solution->getNbValues() ; i++){
                solution.push_back(face_in_solution->value(i));
            }
            std::string s_string = solution_to_string(solution);
            return solution;
        }

        TCellID select_Candidate(std::set<std::pair<TCellID,TCellID>> cdt,std::vector<float> pheromones, std::vector<std::pair<TCellID, float>> &face_proba) {
            assert(cdt.size() != 0);
            std::vector<std::pair<TCellID,TCellID>> candidates(cdt.begin(),cdt.end());
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
                if (same_orientation(c.first, c.second)) {
                    heuristic.push_back(1.0);
                } else {
                    heuristic.push_back(0.5);
                }
            }
            for (std::size_t i = 0; i < candidates.size(); i++) {
                sum_all += pow(pheromones[candidates[i].first], alpha) * pow(heuristic[i], beta);
            }

            for (auto h: heuristic) {
                heuristic_sum += pow(h, beta);
            }
            std::cout << "heuristic_sum " << heuristic_sum << " " << alpha << " " << beta << std::endl;
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
            assert( selected.first < m.getNbFaces());
            return selected.first;
        }

        bool same_orientation(TCellID face_id, TCellID edge_id){
            bool same_orientation = false;
            gmds::Edge edg = m.get<gmds::Edge>(edge_id);
            std::vector<gmds::TCellID> faces;
            edg.getIDs<gmds::Face>(faces);
            TCellID face_solution = 0;// bad
            int count = 0;
            // Where are our marked faces
            for (auto f: faces) {
                if (this->face_in_solution->value(f) != 0) {
                    face_solution = f;
                    break;
                }
            }
            //assert(face_id > solution.size());
            gmds::Face f1 = m.get<gmds::Face>(face_solution), f2 = m.get<gmds::Face>(face_id);
            std::cout << "#FACES " << f1.id() << " " << f2.id() << std::endl;
            if (f1.normal() == f2.normal() || f1.normal() == f2.normal().opp()) {
                same_orientation = true;
            }
            return same_orientation;

        }
        /** Return true if f1 f2 have a region in common.
         * */
        bool edge_turn_topo(TCellID f1, TCellID f2){
            Face face1,face2;
            std::vector<TCellID> r1,r2,res;
            face1 = m.get<Face>(f1);
            face2 = m.get<Face>(f2);
            r1 = face1.getIDs<Region>();
            r2 = face2.getIDs<Region>();
            std::sort(r1.begin(),r1.end());
            std::sort(r2.begin(),r2.end());
            std::set_intersection(r1.begin(), r1.end(),r2.begin(),r2.end(), res.begin());
            if( res.empty() ){
                return false;
            }else{
                return true;
            }
        }

        /**
         * Return true if the
        */
        bool is_vertex_valid( TCellID node ){
            bool res;
            Node n = m.get<Node>(node);
            std::vector<TCellID> edges = n.getIDs<Edge>(),face_in_s;
            for(auto e : edges){
                Edge edge = m.get<Edge>(e);
                std::vector<TCellID> faces = edge.getIDs<Face>(),in_solution;
                for(auto f : faces){
                    if(face_in_solution->value(f) != 0){
                        in_solution.push_back(f);
                    }
                }
                // Différents cas selon position du sommet
                //nb_faces_in_solutions(e)
            }
            return res;
        }

        /* Return 1 if it's ok to have this face as a candidate and 0 if it's not
         * 2 faces in solution for the edge already
         * */
        int edge_of_faces_solution(TCellID face){
            gmds::Face f = m.get<gmds::Face>(face);
            std::vector<TCellID> edges;
            f.getIDs<gmds::Edge>(edges);
            for(auto e : edges){
                if(faces_in_solutions(e) > 1){
                    return 0;
                }
            }
            return 1;
        }

        int faces_in_solutions(TCellID edge) {
            gmds::Edge edges = m.get<gmds::Edge>(edge);
            std::vector<TCellID> faces = edges.getIDs<gmds::Face>();
            int cpt = 0;
            for (auto f: faces) {
                if (this->face_in_solution->value(f)) {
                    cpt++;
                    std::cout << f << std::endl;
                }
            }
            return cpt;
        }

        /** **/

        // Faces constraints
        // Bloc existance
        Env(TInt x_n, TInt y_n, TInt z_n, const int *exist_tens);

        void initialize_hex_number_edge() {
            for (auto e: m.edges()) {
                Edge edge = m.get<Edge>(e);
                this->hex_number_edge->set(e, edge.nbRegions());
            }
        }

        void initialize_hex_number_face() {
            for (auto f: m.faces()) {
                Face face = m.get<Face>(f);
                std::vector<TCellID> regions;
                face.getIDs<Region>(regions);
                this->hex_number_face->set(f, regions.size());
            }
        }

        void do_classification_edge();

        void do_classification_face();

        void init_disable_faces(){
            disable_faces = std::vector<char>(m.getNbFaces(),0);
        }

        /* Disable faces of a candidate
         * This method is bad we could do that in a better way using only one structure
         * */
        void disable_faces_of_selected(std::vector<std::pair<gmds::TCellID, gmds::TCellID>> candidates, std::pair<TCellID,TCellID> selected){
            TCellID edge = selected.second;
            for(auto x : candidates ){
                if(x.second == edge){
                    disable_faces[x.first] = 1;
                }
            }
        }

        /** Read a vtk file and set if a face is in the solution or not
         * Still not working
         * **/
        Env(std::string filename);

        /** alpha,
         *  beta,
         *  p,
         *  tmin,
         *  tmax,
         *  nb_ants
         * */
        void setup_intial_front_and_edges(std::string filename) {
            // File pointer
            std::fstream fin;
            std::string line;
            std::vector<std::string> DataStore;

            std::ifstream file(filename);

            while (getline(file, line, ',')) {
                DataStore.push_back(line);
                getline(file, line);
            }

            for (unsigned int i = 0; i < DataStore.size(); i++) {
                std::cout << DataStore[i] << std::endl;
                this->initial_faces.push_back(stoi(DataStore[i]));
            }
            std::cout << "#initial_face" << initial_faces.size() << std::endl;
            std::cout << "#initial_front" << initial_front.size() << std::endl;

            this->initial_front = constructInitialEdgeFront_v2(this->initial_faces);
            assert(!initial_faces.empty() && !initial_front.empty());
            //this->initial_front = constructInitialEdgeFront_v2(this->initial_faces);
        }

        void setup_intial_front_and_edges(std::vector<std::pair<TCellID, TCellID>> faces) {
            this->initial_front = constructInitialEdgeFront(faces);
            std::set<TCellID> fac;
            std::set<gmds::TCellID> ini_faces;
            for (auto s: faces) {
                ini_faces.insert(s.first);
                ini_faces.insert(s.second);
            }
            this->initial_faces = std::vector<TCellID>(ini_faces.begin(), ini_faces.end());
        }


        void
        setup_ant_colony_params(int alpha, int beta, float p, float tmin, float tmax, int nb_ants, int iter_global,int n_run=1,int nb_elite=1,float tau_initial=1.0);

        ~Env(){}

        std::vector<std::pair<std::vector<char>,float>> get_k_best_solutions(unsigned int k,std::vector<std::vector<char>> solutions, std::vector<float> quality){
            std::cout << "get_k_best_solutions" << std::endl;
            std::vector<std::pair<std::vector<char>,float>> k_best;
            // making pair
            std::vector<std::pair<std::string, float>> order(solutions.size());
            for (auto it = solutions.begin(); it != solutions.end(); ++it) {
                order[std::distance(solutions.begin(),it)] = std::make_pair(solution_to_string(*it), quality[std::distance(solutions.begin(),it)]);
            }
            std::cout << "Sorting" << std::endl;
            std::sort(order.begin(), order.end(), ordering_greater());
            for (auto it: order){
                std::cout << it.first << " (" << it.second << ") "<< std::endl;
            }
            std::cout << std::endl;
            order.erase( unique( order.begin(), order.end() ), order.end() );
            for (auto it: order){
                std::cout << it.first << " (" << it.second << ") " << std::endl;

            }
            std::cout << std::endl;
            // Converting back to vector...
            for(std::size_t i = 0; ( i < k ) && (i < order.size() ) ; i++){
                k_best.push_back(std::make_pair(std::vector<char>(order[i].first.begin(), order[i].first.end()),order[i].second));
                std::transform(std::begin(k_best[i].first),std::end(k_best[i].first),std::begin(k_best[i].first),[](int x){return x-'0';});
            }
            std::cout << "Solution" << std::endl;
            for(auto x : k_best){
                for(auto y : x.first){
                    std::cout << y;
                }
                std::cout << std::endl;
            }
            std::cout << "Solutoins size " << solutions.size() << " Unique solutions " << k_best.size() << std::endl;
            return k_best;
        }

        std::vector<char> string_to_solution(std::string sol ){
            std::vector<char> v(sol.begin(),sol.end());
            for(auto & it: v){
                char c = it;
                it = c;
            }
            return v;
        }

        void setup_solution_from_string(std::string solution) {
            std::vector<char> v_solution;

            for (std::size_t i = 0; i < solution.size(); i++) {
                // Simply convert char '0' or '1' to int
                v_solution.push_back(static_cast<int>(solution[i]) - '0');
            }
            std::cout << v_solution.size() << " " << m.getNbFaces() << std::endl;
            assert(v_solution.size() == m.getNbFaces());
            std::vector<TCellID> f_ini;
            writeSolution(v_solution, "default.vtk", {},{},{});

            char *command = (char *) "paraview --script=pv.py";
            system(command);
        }

        /* Write vtk format fname [*.vtk] */
        void writeVTK(std::string fname);

        /*Write solution in VTK*/
        void writeSolution(std::vector<char> solution, std::string fname,
                           std::vector<float> pheromones, std::vector<TCellID> history_faces, std::vector<TCellID> history_edges);


        std::vector<TCellID> solution_to_faces(std::vector<char> solution) {
            std::vector<TCellID> faces;
            for (size_t i = 0; i < solution.size() - 1; i++) {
                if (solution[i] == 1) {
                    faces.push_back(i);
                }
            }
            return faces;
        }

        void setVerbosity(int verbosity) {
            this->verbosity = verbosity;
        }

        /** Return Edges of a front i.e vector of faces*/
        std::vector<gmds::TCellID> getEdgesOfFront(std::vector<gmds::TCellID> front, std::vector<char> solution);

        /** Return Candidates faces of a front i.e vector of edges*/
        std::vector<gmds::TCellID> getFacesCandidates(std::vector<gmds::TCellID> edges, std::vector<char> solution);

        std::vector<double> heuristicFactor(std::vector<gmds::TCellID> candidates, std::vector<char> solution) {
            std::vector<double> heuristic;
            std::vector<int> sh;
            int shared_with_front;
            sh.reserve(candidates.size());
            heuristic.reserve(candidates.size());
            std::cout << "shared with front " << std::endl;
            for (auto c: candidates) {
                shared_with_front = sharedEdgeWithFront(c, solution);
                std::cout << shared_with_front << ' ';
                sh.push_back(shared_with_front);
            }
            float sum = 0;
            for (auto &it: sh) {
                sum += it;
            }
            std::vector<std::pair<float, float>> proba;
            proba.reserve(candidates.size());
            float prev = 0.0f, cur;
            for (auto it = candidates.begin(); it != candidates.end(); ++it) {
                heuristic.push_back(float(sh[std::distance(candidates.begin(), it)]) / sum);
            }
            /*
            for(auto it = candidates.begin(); it != candidates.end() ; ++it){
                cur = prev + float(sh[std::distance(candidates.begin(),it)])/sum;
                proba.push_back(std::pair<float,float>(prev,cur));
                prev = cur;
            }
            */
            return heuristic;
        }

        void status() {
            std::cout << "#Nodes " << m.getNbNodes() << std::endl;
            std::cout << "#Edges " << m.getNbEdges() << std::endl;
            std::cout << "#Faces " << m.getNbFaces() << std::endl;
            std::cout << "#Regions " << m.getNbRegions() << std::endl;
        }

        /** Return number of edge of 'face_id' in common with the front*/
        int sharedEdgeWithFront(TCellID face, std::vector<char> solution) {
            int cpt = 0;
            gmds::Face f = m.get<gmds::Face>(face);
            std::vector<Edge> edges;
            f.get<gmds::Edge>(edges);
            std::vector<Face> fb;
            //std::cout << "#edge size " <<  edges.size() << std::endl;
            for (auto e: edges) {
                e.get<gmds::Face>(fb);
                for (auto f_id: fb) {
                    //if(( this->front->value( f_id.id() ) > FaceStatus::Candidate ) && ( f_id.id() != face )){
                    if ((solution[f_id.id()] == 1) && (f_id.id() != face)) {
                        cpt++;
                    }
                }
            }
            //std::cout << "#cpt " << cpt << std::endl;
            return cpt;
        }

        /** Return a candidate from a vector*/

        std::pair<gmds::TCellID, gmds::TCellID>
        selectCandidate(std::vector<std::pair<gmds::TCellID, gmds::TCellID>> candidates, std::vector<char> solution,
                            std::vector<float> pheromones, std::vector<std::pair<TCellID, float>> &proba);

        /** Update front */
        void updateFront(std::vector<gmds::TCellID> front);

        void updateHeuristique(std::vector<double> heuristic, std::vector<TCellID> candidates) {
            this->heuristic->value(0.0f);
            for (auto it = heuristic.begin(); it != heuristic.end(); ++it) {
                this->heuristic->set(candidates[std::distance(heuristic.begin(), it)], *it);
            }

        }

        void updatePheromones(std::vector<float> pheromone) {
            for (auto it = pheromone.begin(); it != pheromone.end(); ++it) {
                this->pheromone->set(std::distance(pheromone.begin(), it), *it);
            }
        }
        //template <class T> void updateVariable(Variable<T> variable, std::vector<T> vec);
        /** Update front*/

        void updateCandidates(std::vector<gmds::TCellID> candidates);

        void write_R(std::string filename) {
            IGMeshIOService ios(&this->m);
            VTKWriter writer(&ios);
            writer.setCellOptions(N | E | R);
            writer.setDataOptions(N | E | R);
            writer.write(
                    filename); // paraview ctrl + space (Shrink) pour avoir uniquement les hex (enlever les F dans l'écriture)
        }

        std::list<gmds::TCellID>
        getFirstEdgeFront(std::vector<gmds::TCellID> front_initial, std::vector<char> solution);

        void update_probabilities(std::vector<std::pair<TCellID, float>> proba) {
            std::cout << "update_proba" << std::endl;

            for (auto f_p: proba) {
                assert(f_p.first < m.getNbFaces());
                std::cout << f_p.first << f_p.second << std::endl;
                this->probabilities->set(f_p.first, f_p.second);
            }
        }


        /** Utils*/
        void show_tcell_vector(std::vector<gmds::TCellID> v);

        void resetVariables();

        void updateFrontInitial(std::vector<TCellID> FrontInitial);

        /*Return Edges opposite of a face TCELLID*/
        std::vector<TCellID> getEdgesOpposite(TCellID face);

        void show_normals(std::vector<TCellID> faces);

        int evaluate(std::vector<char> solution) {
            resetVariables();

            for (gmds::FaceContainer::iterator it = m.faces().begin(); it != m.faces().end(); ++it) {
                /**if (solution[std::distance(m.faces.begin(), it)] == 1) {
                    front->set(*it, FaceStatus::Front);
                }
                */
            }
            for (auto it = solution.begin(); it != solution.end(); ++it) {
                front->set(*it, FaceStatus::Initial);
            }
            return 1;//numberEdgeTurn(solution);

        }

        int numberEdgeTurn(std::vector<char> solution) {
            int counter = 0;
            for (auto e_id: m.edges()) {
                if (edgeTurn(e_id, solution))
                    counter++;
            }
            return counter;
        }

        /*!
         * Execute algorithm
         * front_initial := front of faces initialy constraint
         * */
        void execute();

        int getNbFaces() {
            return m.getNbFaces();
        }

        bool edgeTurn(TCellID edge, std::vector<char> solution) {
            bool edge_turn = true;
            gmds::Edge edg = m.get<gmds::Edge>(edge);
            std::vector<gmds::TCellID> faces;
            edg.getIDs<gmds::Face>(faces);
            TCellID fac[2];
            int count = 0;
            // Where are our marked faces
            for (auto f: faces) {
                //if (this->front->value(f) > FaceStatus::Candidate) {
                if (solution[f] != 0) {
                    fac[count] = f;
                    count++;
                }
            }
            // assert(count < 3);
            if( count == 2 && solution[fac[0]] == 3 && solution[fac[1]] == 3 ){
                return false;
            }else if (count == 2 ) { // 2 faces in front
                gmds::Face f1 = m.get<gmds::Face>(fac[0]), f2 = m.get<gmds::Face>(fac[1]);
                if (f1.normal() == f2.normal() || f1.normal() == f2.normal().opp()) {
                    edge_turn = false;
                }
            } else {
                edge_turn = false;
            }
            return edge_turn;
        }

        /*Return true if a 'face_id' have the same normal (or opposite normal) of a face solution contained in solution
         **/
        bool same_orientation(TCellID face_id, TCellID edge_id, std::vector<char> solution) {
            bool same_orientation = false;
            gmds::Edge edg = m.get<gmds::Edge>(edge_id);
            std::vector<gmds::TCellID> faces;
            edg.getIDs<gmds::Face>(faces);
            TCellID face_solution = 0;// bad
            int count = 0;
            // Where are our marked faces
            for (auto f: faces) {
                if (solution[f] == 1) {
                    face_solution = f;
                    break;
                }
            }

            //assert(face_id > solution.size());
            gmds::Face f1 = m.get<gmds::Face>(face_solution), f2 = m.get<gmds::Face>(face_id);
            std::cout << "#FACES " << f1.id() << " " << f2.id() << std::endl;
            if (f1.normal() == f2.normal() || f1.normal() == f2.normal().opp()) {
                same_orientation = true;
            }
            return same_orientation;
        }


        void test_edge_to_add() {
            for (auto f: m.faces()) {
                Face face = m.get<Face>(f);
                std::vector<TCellID> edges;
                face.getIDs<Edge>(edges);
                for (auto e: edges) {

                }
            }
        }

        edge_classification get_edge_classification(TCellID edge) {
            return edge_classification::convex;
        }

        face_classification get_face_classification(TCellID face) {
            return face_classification::surface;
        }

        void edges_to_add() {
            std::cout << m.getNbFaces() << std::endl;
            for (auto f: m.faces()) {
                Face face = m.get<Face>(f);
                std::vector<TCellID> edges;
                face.getIDs<Edge>(edges);
                std::cout << edges.size() << std::endl;
                for (auto e: edges) {
                    Edge edg = m.get<Edge>(e);
                    std::vector<Region> reg;
                    std::vector<Region> reg_face;
                    face.get<Region>(reg_face);
                    edg.get<Region>(reg);
                    std::cout << f << "," << e << ": " << classification_face->value(f) << ", "
                              << classification_edge->value(e) << " ( "
                              << edge_to_add(static_cast<face_classification>(classification_face->value(f)),
                                             static_cast<edge_classification>(classification_edge->value(e))) << ")"
                              << std::endl;
                }
            }
        }


        /** Do we have to add 'e' to the front knowing face 'f' ?
         *
         * */
        bool edge_to_add(face_classification f, edge_classification e) {
            bool b;
            if (f == face_classification::surface && e == edge_classification::concave) {
                b = true;
            } else if (f == face_classification::surface && e == edge_classification::convex) {
                b = false;
            } else if (f == face_classification::surface && e == edge_classification::surface) {
                b = true;
            } else if (f == face_classification::surface && e == edge_classification::volume) {
                // impossible
                assert(f == face_classification::surface && e == edge_classification::volume);
            } else if (f == face_classification::surface && e == edge_classification::volume) {
                b = true;
            } else if (f == face_classification::volume && e == edge_classification::concave) {
                b = true;
            } else if (f == face_classification::volume && e == edge_classification::convex) {
                //impossible
                std::cout << "error " << int(f) << "," << int(e) << std::endl;
                assert(1 == 2);
                b = false;
            } else if (f == face_classification::volume && e == edge_classification::surface) {
                b = false;
            } else if (f == face_classification::volume && e == edge_classification::volume) {
                b = true;
            } else {
                std::cout << "error " << int(f) << "," << int(e) << std::endl;

                assert(1 == 2);
            }
            return b;
        }

        /***    Return Edges of a face and his opposite
        *       not processed
        *       face in solution == 1
        ***/
        std::vector<std::pair<TCellID, TCellID>>
        getEdges_v2(TCellID face_id, std::vector<char> processed, std::vector<char> solution) {
            gmds::Face face = m.get<gmds::Face>(face_id);
            std::vector<gmds::Edge> edges;
            face.getOrderedEdges(edges);// get opposite edges
            std::vector<std::pair<TCellID, TCellID>> res, opposites;    // Edge,opposite
            opposites.push_back(std::pair<TCellID, TCellID>(edges[0].id(), edges[2].id()));
            opposites.push_back(std::pair<TCellID, TCellID>(edges[1].id(), edges[3].id()));
            opposites.push_back(std::pair<TCellID, TCellID>(edges[2].id(), edges[0].id()));
            opposites.push_back(std::pair<TCellID, TCellID>(edges[3].id(), edges[1].id()));
            for (auto e = edges.begin(); e != edges.end(); ++e) {
                int count = 0;
                std::vector<TCellID> incident_faces;
                (*e).getIDs<gmds::Face>(incident_faces);
                for (auto face: incident_faces) {
                    if (solution[face] == 1) {
                        count++;
                    }
                }
                //! Condition here: not processed and only one face in solution
                if (processed[(*e).id()] == false && count == 1) {
                    res.push_back(opposites[std::distance(edges.begin(), e)]);
                }
            }
            return res;
        }

        std::vector<std::pair<TCellID,TCellID>> getEdges_v3(TCellID face_id,TCellID edge_id,std::vector<char> processed_edge){
            std::vector<std::pair<TCellID,TCellID>> res;
            Face face = m.get<Face>(face_id);
            std::vector<TCellID> edges;
            face.getIDs<Edge>(edges);
            for( auto e : edges){
                if( processed_edge[e] == 0 && edge_id != e && edge_to_add(static_cast<face_classification>(classification_face->value(face_id)),static_cast<edge_classification>(classification_edge->value(e)))){
                    res.push_back(std::make_pair(e,face_id));
                }
            }
            return res;
        }

        int faces_in_solutions(TCellID edge, std::vector<char> solution) {
            gmds::Edge edges = m.get<gmds::Edge>(edge);
            std::vector<TCellID> faces = edges.getIDs<gmds::Face>();
            int cpt = 0;
            for (auto f: faces) {
                if (solution[f] != 0 ) {
                    cpt++;
                    std::cout << f << std::endl;
                }
            }
            return cpt;
        }

        TCellID face_id_in_solutions(TCellID edge, std::vector<char> Solution) {
            gmds::Edge edges = m.get<gmds::Edge>(edge);
            std::vector<TCellID> faces = edges.getIDs<gmds::Face>();
            TCellID face;
            bool found = false;
            for (auto f: faces) {
                if (Solution[f] != 0) {
                    face = f;
                    found = true;
                }
            }
            assert(found == true);
            return face;
        }

        std::vector<std::pair<TCellID, TCellID>> constructInitialEdgeFront_v2(std::vector<TCellID> initial_faces) {
            std::vector<std::pair<TCellID, TCellID>> initial_front;
            std::vector<std::vector<char>> adj(initial_faces.size(), std::vector<char>(initial_faces.size()));
            for (std::size_t i = 0; i < initial_faces.size(); i++) {
                for (std::size_t j = 0; j < initial_faces.size(); j++) {
                    gmds::Face f1 = m.get<Face>(initial_faces[i]), f2 = m.get<Face>(initial_faces[j]);
                    std::vector<TCellID> nodes = m.getCommonNodes(f1, f2);
                    if ((i != j) && adj[i][j] == 0 && !nodes.empty()) {
                        TCellID initial_edge = getCommonEdge(initial_faces[i], initial_faces[j]);
                        std::vector<gmds::Edge> edges_f1, edges_f2;
                        f1.getOrderedEdges(edges_f1);
                        f2.getOrderedEdges(edges_f2);
                        auto e0 = std::find(edges_f1.begin(), edges_f1.end(), m.get<gmds::Edge>(initial_edge));
                        Edge opposite_f1, opposite_f2; // Opposite is  ( x+1 ) mod size ) +2 -1
                        opposite_f1 = edges_f1[((std::distance(edges_f1.begin(), e0) + 3) % edges_f1.size()) - 1];
                        e0 = std::find(edges_f2.begin(), edges_f2.end(), m.get<gmds::Edge>(initial_edge));
                        opposite_f2 = edges_f2[((std::distance(edges_f2.begin(), e0) + 3) % edges_f2.size()) - 1];
                        assert(opposite_f1.id() != opposite_f2.id());
                        initial_front.push_back(std::pair<TCellID, TCellID>(opposite_f1.id(), initial_edge));
                        initial_front.push_back(std::pair<TCellID, TCellID>(opposite_f2.id(), initial_edge));
                    }
                }
            }
            assert(!initial_front.empty());
            return initial_front;
        }

        /** Get if a face have an edge with faces in solution 0 if ok to be a candidate 1 else*/
        int edge_of_faces_solution(TCellID face,std::vector<char> solution){
            gmds::Face f = m.get<gmds::Face>(face);
            std::vector<TCellID> edges;
            f.getIDs<gmds::Edge>(edges);
            for(auto e : edges){
                if(faces_in_solutions(e,solution) > 1){
                    return 0;
                }
            }
            return 1;
        }

        std::pair<TCellID,TCellID> random_selection(std::map<TCellID,std::vector<TCellID>> candidates){
            std::vector<TCellID> v;
            v.reserve(candidates.size());
            for(auto c : candidates){
                v.push_back(c.first);
            }
            int min = 0, max = candidates.size() - 1;
            int randNum = rand() % (max - min + 1) + min;
            return std::make_pair(v[randNum],candidates[v[randNum]][0]);
        }

        /** Return edge id in common with face f1 and f2
         * */
        TCellID getCommonEdge(TCellID f1, TCellID f2) {
            std::vector<TCellID> nodes = m.getCommonNodes(m.get<gmds::Face>(f1), m.get<gmds::Face>(f2));
            gmds::Node n = m.get<gmds::Node>(nodes[0]);
            std::vector<gmds::Edge> edges;
            n.get<gmds::Edge>(edges);
            TCellID target_edge;
            for (auto e: edges) {
                std::vector<gmds::TCellID> no = e.getIDs<gmds::Node>();
                if (no[0] == nodes[1] || no[1] == nodes[1]) {
                    target_edge = e.id();
                    break;
                }
            }
            return target_edge;
        }

        // We use pair here edge, and his opposite
        std::vector<std::pair<TCellID, TCellID>>
        constructInitialEdgeFront(std::vector<std::pair<TCellID, TCellID>> initial_faces) {
            std::vector<std::pair<TCellID, TCellID>> initial_front;
            for (auto pair: initial_faces) {
                TCellID initial_edge = getCommonEdge(pair.first, pair.second);
                // Get edges from f1 f2
                gmds::Face f1 = m.get<gmds::Face>(pair.first), f2 = m.get<gmds::Face>(pair.second);
                std::vector<gmds::Edge> edges_f1, edges_f2;
                f1.getOrderedEdges(edges_f1);
                f2.getOrderedEdges(edges_f2);

                auto e0 = std::find(edges_f1.begin(), edges_f1.end(), m.get<gmds::Edge>(initial_edge));
                Edge opposite_f1, opposite_f2; // Opposite is  ( x+1 ) mod size ) +2 -1
                opposite_f1 = edges_f1[((std::distance(edges_f1.begin(), e0) + 3) % edges_f1.size()) - 1];
                e0 = std::find(edges_f2.begin(), edges_f2.end(), m.get<gmds::Edge>(initial_edge));
                opposite_f2 = edges_f2[((std::distance(edges_f2.begin(), e0) + 3) % edges_f2.size()) - 1];

                assert(opposite_f1.id() != opposite_f2.id());

                initial_front.push_back(std::pair<TCellID, TCellID>(opposite_f1.id(), initial_edge));
                initial_front.push_back(std::pair<TCellID, TCellID>(opposite_f2.id(), initial_edge));
            }
            return initial_front;
        }

        std::vector<char> constructInitialSolution(std::vector<TCellID> initial_faces) {
            std::vector<char> solution(m.getNbFaces());
            std::fill(solution.begin(), solution.end(), 0);
            for (auto f: initial_faces) {
                solution[f] = 1;
            }
            return solution;
        }

        std::string solution_to_string(std::vector<char> solution) {
            std::string s;
            for (auto x: solution) {
                s += std::to_string(x);
            }
            return s;
        }

        std::vector<char> run_ant(std::vector<float> pheromones);
        //! get opposite edge verifier et adapter pour obtenir l'arrete opposer dans une face a partir d'une arrete
        /** Writer
         * - Write in dir ...
         * - Dir 0...
         *      0-end_edge_turn.vtk
         *          pheromones: float
         *          heuristique: float
         *          probabilities: float
         *          Face status: Solution, Candidate, Nothing
         *      .txt with face selected
         *          For each turn
         *              edge in front, candidates for this edges,
         *              candidates array
         *              candidates boolean
         *              pheromon array
         *              heuristic array
         *              candidate selected
         *              edge removed
         *              edges inserted
         **/

        bool is_selection_valid(std::vector<char> solution) {
            bool is_valid = true;
            std::cout << "Nb cells: " << m.getNbRegions() << std::endl;
            std::cout << "Nb faces: " << m.getNbFaces() << std::endl;
            std::cout << "Nb edges: " << m.getNbEdges() << std::endl;
            std::cout << "Nb nodes: " << m.getNbNodes() << std::endl;
            std::vector<TCellID> E1H, E2H, E3H, E4H;

            for (auto edge_id: m.edges()) {
                Edge e = m.get<Edge>(edge_id);

                std::vector<TCellID> regions_of_e = e.getIDs<Region>();
                if (regions_of_e.size() == 1) {
                    E1H.push_back(edge_id);
                } else if (regions_of_e.size() == 2) {
                    E2H.push_back(edge_id);
                } else if (regions_of_e.size() == 3) {
                    E3H.push_back(edge_id);
                } else {
                    E4H.push_back(edge_id);
                }
            }
            std::cout << "|E1H|=" << E1H.size() << std::endl;
            std::cout << "|E2H|=" << E2H.size() << std::endl;
            std::cout << "|E3H|=" << E3H.size() << std::endl;
            std::cout << "|E3H|=" << E4H.size() << std::endl;
            assert((E1H.size() + E2H.size() + E3H.size() + E4H.size()) == m.getNbEdges());

            for (auto e: E4H) {
                gmds::Edge edg = m.get<Edge>(e);
                if (faces_in_solutions(e, solution) % 2 != 0) {
                    is_valid = false;
                    break;
                }
            }


            return is_valid;
        }

        float best_quality(std::vector<char> solution) {
            TInt nb_node = m.getNbNodes();
            TInt nb_edge = m.getNbEdges();
            TInt nb_face = m.getNbFaces();
            return float(nb_face + nb_edge + nb_node);
        }

        float solution_quality(std::vector<char> solution) {
            TInt nb_node = m.getNbNodes();
            TInt nb_edge = m.getNbEdges();
            TInt nb_face = m.getNbFaces();
            int total_edges = 0;
            for (auto e: m.edges()) {
                if (faces_in_solutions(e, solution) == 2)
                    total_edges += 1;
            }
            std::set<TCellID> vertex_set;
            for(auto f : m.faces() ) {
                if (solution[f] != 0){
                    gmds::Face face = m.get<Face>(f);
                    std::vector<TCellID> nodes;
                    face.getIDs<Node>(nodes);
                    for(auto n : nodes){
                        vertex_set.insert(n);
                    }
                }
            }


            float edge_part = static_cast< float >(total_edges - nb_edge_turn(solution)) /
                              static_cast< float >(total_edges);
            float vertex_part = static_cast<float >(nb_node - nb_vertex_turn(solution)) /static_cast<float>(nb_node) ;
            std::size_t nb_element = std::count(solution.begin(), solution.end(), 1);

            float element_part = static_cast<float >(nb_face - nb_element) / static_cast<float >(nb_face);

            if(vertex_part < 0.1) {
                vertex_part = 0.1;
            }
            //float res=edge_part + vertex_part + element_part;
            std::cout << "Recap" << std::endl;
            std::cout << total_edges << " " << nb_edge_turn(solution) << std::endl ;
            std::cout << nb_node << " " << nb_vertex_turn(solution) << std::endl ;
            std::cout << edge_part << " " << element_part << std::endl;
            float res=30*(edge_part+0.1*element_part);
            return res;
        }

        bool is_solution_valid(std::vector<char> solution) {
            std::vector<TCellID> edge_border;
            for (auto e: m.edges()) {
                if (faces_in_solutions(e, solution) == 1) {
                    edge_border.push_back(e);
                }
            }
            for (auto e: edge_border) {
                TCellID face = face_id_in_solutions(e, solution);
                gmds::Edge edge = m.get<gmds::Edge>(e);

                if (this->hex_number_edge->value(e) > 3) {
                    return false;
                }
                if ((face_boundary->value(face) == 1) && (this->hex_number_edge->value(e) != 1)) {
                    return false;
                }
            }
            return true;
        }

        int nb_vertex_turn(std::vector<char> solution) {
            std::vector<char> vtx_turn;
            int vtx_turn_count = 0;
            for (auto n: m.nodes()) {
                bool temp = is_vertex_turn(n, solution);
                if (temp) {
                    vtx_turn.push_back(temp);
                    vtx_turn_count++;
                }
            }
            return vtx_turn_count;
        }

        int nb_edge_turn(std::vector<char> solution) {
            std::vector<char> edge_turn;
            int edge_turn_count = 0;
            for (auto e: m.edges()) {
                bool temp = edgeTurn(e, solution);
                if (temp) {
                    edge_turn.push_back(temp);
                    edge_turn_count++;
                }
            }
            return edge_turn_count;
        }

        bool is_vertex_turn(TCellID node_id, std::vector<char> solution){
            std::cout << node_id << std::endl;
            bool res = false;
            gmds::Node n = m.get<Node>(node_id);
            std::vector<TCellID> edges = n.getIDs<Edge>();
            std::vector<TCellID> nodes;
            std::vector<TCellID> nodes_tmp;
            std::vector<math::Point> points_a;
            std::vector<math::Point> points_b;
            std::vector<math::VectorDyn> vdyn;
            for (auto e: edges) {
                Edge edge = m.get<Edge>(e);
                nodes_tmp = edge.getIDs<Node>();
                Node n0 = m.get<Node>(nodes_tmp[0]);
                Node n1 = m.get<Node>(nodes_tmp[1]);
                points_a.push_back(n0.point());
                points_b.push_back(n1.point());
                math::VectorDyn a = math::VectorDyn(n0.X() - n1.X(), n0.Y() - n1.Y(), n0.Z() - n1.Z());
                vdyn.push_back(a);
            }
            int i = 0;
            while (!res && i < vdyn.size()) {
                for (int j = i; j < vdyn.size() - 1; j++) {
                    // vectors have to be colinear and the value of edge turn different
                    if (vdyn[i].isColinear(vdyn[j]) && (edgeTurn(edges[i], solution) ^ edgeTurn(edges[j],
                                                                                                solution))) {//&& (is_edge_turn(edges[i]) ^ is_edge_turn(edges[j])
                        std::cout << " is_edg_turn " << edgeTurn(edges[i], solution) << edgeTurn(edges[j], solution)
                                  << std::endl;
                        res = true;
                    }
                }
                i++;
            }
            return res;
        }
    };
}
#endif //ENV_H
