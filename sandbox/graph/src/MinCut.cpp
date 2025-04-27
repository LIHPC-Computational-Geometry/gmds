//
// Created by calderans on 23/07/20.
//

#include <iostream>

#include "gmds/graph/MinCut.h"
#include "../../external/gco-v3/inc/GCoptimization.h"

using namespace gmds;
using namespace graph;

auto comp = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {return a.second > b.second; };

MinCut::MinCut(Mesh* AMesh):m_mesh(AMesh) {
}
/*----------------------------------------------------------------------------*/

void MinCut::graph_cut(Variable<int> *ATetAssign, Variable<int> *ATetResult,
                       Variable<double> *AFaceWeight) {

    int nb_tets = m_mesh->getNbRegions();
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(nb_tets,2);

    for(auto fid:m_mesh->faces()){
        Face f = m_mesh->get<Face>(fid);
        std::vector<TCellID > adj_r = f.getIDs<Region>();
        if(adj_r.size()==2){
            math::Vector3d n = f.normal();
            gc->setNeighbors(adj_r[0],adj_r[1], AFaceWeight->value(f.id()));
        }

    }
    // first set up the array for data costs
    double *data = new double[nb_tets*2]; //2 values 0 (left)/(1 right)

    for(auto rid:m_mesh->regions()){

        if(ATetAssign->value(rid)==1){
            data[rid * 2] = 0;
            data[rid * 2 + 1] = 1000;

        }
        else if (ATetAssign->value(rid)==2){
            data[rid * 2] = 1000;
            data[rid * 2 + 1] = 0;
        }
        else {
            //general case for inner tet
            data[rid * 2] = 1;
            data[rid * 2 + 1] = 1;
        }
    }
    gc->setDataCost(data);

    int *result = new int[nb_tets];   // stores result of optimization

    // next set up the array for smooth costs
    double *smooth = new double[2*2];
    for ( int l1 = 0; l1 < 2; l1++ ) {
        for (int l2 = 0; l2 < 2; l2++) {
            if (l1 == l2) {
                smooth[l1 + l2 * 2] = 0;
            } else {
                smooth[l1 + l2 * 2] = 100;
            }
        }
    }
    gc->setSmoothCost(smooth);

	 std::cout<<"Before optimization energy is "<<gc->compute_energy()<<std::endl;
    gc->swap(100);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
	 std::cout<<"After optimization energy is  "<<gc->compute_energy()<<std::endl;
    std::cout<<std::endl;

    for(auto rid:m_mesh->regions()) {
        int label = gc->whatLabel(rid);
        ATetResult->set(rid,label+1);
    }
    delete [] smooth;
    delete [] data;
    delete  gc;

}
/*----------------------------------------------------------------------------*/
void MinCut::graph_cut(std::map<TCellID, int> ATetList, std::vector<TCellID> AFaceList, Variable<int> *ATetAssign,
                       Variable<int> *ATetResult, Variable<double> *AFaceWeight) {
    int nb_tets = ATetList.size();
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(nb_tets,2);

    for(auto fid:AFaceList){
        Face f = m_mesh->get<Face>(fid);
        std::vector<TCellID > adj_r = f.getIDs<Region>();
        if(adj_r.size()==2){
            gc->setNeighbors(ATetList[adj_r[0]],ATetList[adj_r[1]], AFaceWeight->value(f.id()));
        }

    }
    // first set up the array for data costs
    double *data = new double[nb_tets*2]; //2 values 0 (left)/(1 right)

    for(auto rid:ATetList){

        if(ATetAssign->value(rid.first)==1){
            data[rid.second * 2] = 0;
            data[rid.second * 2 + 1] = 10000;

        }
        else if (ATetAssign->value(rid.first)==2){
            data[rid.second * 2] = 10000;
            data[rid.second * 2 + 1] = 0;
        }
        else {
            //general case for inner tet
            data[rid.second * 2] = 1;
            data[rid.second * 2 + 1] = 1;
        }
    }
    gc->setDataCost(data);

    int *result = new int[nb_tets];   // stores result of optimization

    // next set up the array for smooth costs
    double *smooth = new double[2*2];
    for ( int l1 = 0; l1 < 2; l1++ ) {
        for (int l2 = 0; l2 < 2; l2++) {
            if (l1 == l2) {
                smooth[l1 + l2 * 2] = 0;
            } else {
                smooth[l1 + l2 * 2] = 100;
            }
        }
    }
    gc->setSmoothCost(smooth);

    printf("\nBefore optimization energy is %d",gc->compute_energy());
    gc->swap(100);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
    printf("\nAfter optimization energy is %d",gc->compute_energy());
    std::cout<<std::endl;

    for(auto rid:ATetList) {
        int label = gc->whatLabel(rid.second);
        ATetResult->set(rid.first,label+1);
    }
    delete [] smooth;
    delete [] data;
    delete  gc;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> MinCut::shortest_path(TCellID ASource, std::vector<TCellID> ANodeGraph,
                                           std::map<TCellID, double> AEdgeWeight, std::vector<TCellID> ATargetNodes) {
    int NbNodes = ANodeGraph.size(), NbEdges = AEdgeWeight.size();

    std::vector<std::vector<std::pair<int, double>>> G(NbNodes);                                    // The graph is a vector with NbNodes nodes.
    // Each node is connected to others nodes via weighted edges. This information is stored in a vector of pair
    for (auto e : AEdgeWeight) {
        TCellID u_id, v_id;
        int u, v;


        Edge edge = m_mesh->get<Edge>(e.first);
        u_id = edge.getIDs<Node>()[0];
        v_id = edge.getIDs<Node>()[1];
        double w = e.second;



        auto itr = std::find(ANodeGraph.begin(),ANodeGraph.end(), u_id);
        u = std::distance(ANodeGraph.begin(),itr);
        itr = std::find(ANodeGraph.begin(),ANodeGraph.end(), v_id);
        v = std::distance(ANodeGraph.begin(),itr);



        G[u].push_back(std::make_pair(v, w));                                          // Each pair contains : first=index of the node connected to u, second=weight/distance/codst of the path from v to u
        G[v].push_back(std::make_pair(u, w));                                          // Comment this line if the graph use directed edges.
        // With undirected edges, create link from v to u and u to v. Both with weight w
    }

    TCellID StartNode;

    for(int i = 0; i<NbNodes; i++){
        if(ANodeGraph[i] == ASource){
            StartNode = i;
        }
    }

    std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());                   // Distances is a vector of NbNodes cells. All cells are initialized with max()
    // Distances[i] is the distance from StartNode to node whose index is i

    Distances[StartNode] = 0;                                                     // Distance to StartNode is initialized to 0

    std::vector<int> Parents(NbNodes, -1);                                             // Parents is a vector of NbNodes cells. All cells are initialized with -1

    // Priority queue storing pairs and using a specific comparator function
    // Because of the comparator we need to specify the 3 parameters
    // The comparator make sure that the closest node is always on top of the queue
    // Each pair is made of : index of the node and the distance to StartNode
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, decltype(comp)> Q(comp);
    Q.push(std::make_pair(StartNode, 0));                                              // Initialize the priority queue with StartNode

    int v0 = Q.top().first;                                                      // get the index of the nearest node
    int w0 = Q.top().second;                                                     // get the weight/cost/distance
    Q.pop();

    if (w0 <= Distances[v0]) {                                                    // Pay attention to this test.
        // It can be removed, however, it avoid duplicated work

        for (const auto& i : G[v0]) {                                              // v is the index of the nearest node
            auto v2 = i.first;                                                      // For each node connected to node v
            auto w2 = i.second;

            if (Distances[v0] + w2 < Distances[v2]) {                                // If distance from StartNode to v2 thru v is shorter then the current distance from StartNode to v2
                Distances[v2] = Distances[v0] + w2;                                    // then update the distance from StartNode to v2 and parent[v2]
                Parents[v2] = v0;                                                      // https://www.youtube.com/watch?v=8Ls1RqHCOPw
                Q.push(std::make_pair(v2, Distances[v2]));

            }
        }
    }

    while (!Q.empty()) {                                                          // Dijkstra
        auto v = Q.top().first;                                                      // get the index of the nearest node
        auto w = Q.top().second;                                                     // get the weight/cost/distance
        Q.pop();

        math::Point n_v_p = m_mesh->get<Node>(ANodeGraph[Parents[v]]).getPoint();
        math::Point n_v  = m_mesh->get<Node>(ANodeGraph[v]).getPoint();

        math::Vector3d parent_vec(n_v_p.X()+n_v.X(),n_v_p.Y()+n_v.Y(),n_v_p.Z()+n_v.Z());
        parent_vec.normalize();

        if (w <= Distances[v]) {                                                    // Pay attention to this test.
            // It can be removed, however, it avoid duplicated work

            for (const auto& i : G[v]) {                                              // v is the index of the nearest node
                auto v2 = i.first;                                                      // For each node connected to node v
                auto w2 = i.second;

                math::Vector3d new_vec;

                //On calcule l'angle entre l'arête précédente et les nouvelles possibles

                //On crée le vecteur entre v et v2 (l'arc à tester)

                math::Point n_v2 = m_mesh->get<Node>(ANodeGraph[v2]).getPoint();

                new_vec.setXYZ(n_v.X()+n_v2.X(),n_v.Y()+n_v2.Y(),n_v.Z()+n_v2.Z());
                new_vec.normalize();

                double dot = (1-fabs(new_vec.dot(parent_vec)));

                if (Distances[v] + w2 + dot < Distances[v2]) {                                // If distance from StartNode to v2 thru v is shorter then the current distance from StartNode to v2
                    Distances[v2] = Distances[v] + w2 + dot;                                    // then update the distance from StartNode to v2 and parent[v2]
                    Parents[v2] = v;                                                      // https://www.youtube.com/watch?v=8Ls1RqHCOPw
                    Q.push(std::make_pair(v2, Distances[v2]));
                }
            }
        }
    }

    /*for (auto i = 0; i != NbNodes; ++i) {                                          // display the results
        std::cout << "\nPath from node " << ASource << " to node " << ANodeGraph[i] << " cost " << Distances[i] << std::endl;

        std::cout << i;
        for (auto p = Parents[i]; p != -1; p = Parents[p])
            std::cout << " <- " << ANodeGraph[p];                                                      // when links are not bi directional the output is accurate when using <- instead of ->
        std::cout << std::endl;                                                               // otherwise it make no difference
    }*/




    std::map<TCellID,double> targets_weight;
    for(auto t : ATargetNodes){
        auto itr = std::find(ANodeGraph.begin(),ANodeGraph.end(), t);
        int index = std::distance(ANodeGraph.begin(),itr);

        //std::cout << "\nPath from node " << ASource << " to node " << ANodeGraph[index] << " cost " << Distances[index] << std::endl;
        targets_weight.emplace(t,Distances[index]);

        //std::cout << ANodeGraph[index];
        //for (auto p = Parents[index]; p != -1; p = Parents[p])
        //std::cout << " <- " << ANodeGraph[p];                                                      // when links are not bi directional the output is accurate when using <- instead of ->
        //std::cout << std::endl;

    }

    double min_weight = 9999;
    TCellID min_id;
    for(auto t_w : targets_weight){
        if(min_weight > t_w.second && t_w.first != ASource){
            min_weight = t_w.second;
            min_id = t_w.first;
        }
    }
    std::vector<TCellID> result;
    auto itr = std::find(ANodeGraph.begin(),ANodeGraph.end(), min_id);
    int index = std::distance(ANodeGraph.begin(),itr);
    result.push_back(min_id);
    for (auto p = Parents[index]; p != -1; p = Parents[p]){
        result.push_back(ANodeGraph[p]);
    }
    std::cout<<"Nearest point from "<<ASource<<" is "<<min_id<<" with cost "<<min_weight<<std::endl;
    return result;


    //getchar();
}
/*----------------------------------------------------------------------------*/


