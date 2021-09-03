//
// Created by calderans on 10/03/20.
//

#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>


#include <fstream>


#include <gmds/graph/TopologyGraph.h>
/*----------------------------------------------------------------------------*/

using namespace gmds;
using namespace graph;
/*----------------------------------------------------------------------------*/
TopologyGraph::TopologyGraph(cad::FACManager* Amanager, Mesh* Amesh):m_manager(Amanager),m_mesh(Amesh){

    try {
        m_surface_color = m_mesh->newVariable<int, GMDS_FACE>("surface_color");
    }catch (GMDSException e) {
        m_surface_color = m_mesh->getVariable<int, GMDS_FACE>("surface_color");
    }
    try {
        m_isomorph = m_mesh->newVariable<int, GMDS_FACE>("isomorph");
    }catch (GMDSException e) {
        m_isomorph = m_mesh->getVariable<int, GMDS_FACE>("isomorph");
    }
}
/*----------------------------------------------------------------------------*/
bool TopologyGraph::buildGraph() {

    m_graph.resize(m_manager->getNbSurfaces());

    std::cout<<"NB surfaces = "<<m_manager->getNbSurfaces()<<std::endl;
    std::cout<<"graph size "<<m_graph.size()<<std::endl;

    for(int i=0; i<m_graph.size(); i++){
        m_graph[i].resize(m_manager->getNbSurfaces());
        for(auto n2:m_graph[i]){
            n2 = 0;
        }
    }

    std::vector<cad::GeomCurve*> curves;
    m_manager->getCurves(curves);
    for(auto c:curves){

        cad::GeomSurface* s0 = c->surfaces()[0];
        int id0 = s0->id()-1;

        cad::GeomSurface* s1 = c->surfaces()[1];
        int id1 = s1->id()-1;

        //if the curve have only one point or zero it is closed we note the arc with 2
        if(c->points().size()<=1){
            m_graph[id0][id1] = 2;
            m_graph[id1][id0] = 2;
        }else{
            m_graph[id0][id1] = 1;
            m_graph[id1][id0] = 1;
        }

    }

    return true;
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::printGraph() {
    std::cout<<"  ";
    for(int i =1; i<=m_graph.size();i++){
        std::cout<<i<<",";
    }
    std::cout<<std::endl;
    int in=1;
    for(auto n: m_graph){
        std::cout<<in++<<":";
        for(auto n2: n){
            std::cout<<n2<<",";
        }
        std::cout<<std::endl;
    }
    std::cout<<"\nv2"<<std::endl;
    std::cout<<"  ";
    for(int i =1; i<=m_graph.size();i++){
        std::cout<<i<<",";
    }
    std::cout<<std::endl;


    for(int i =0; i<m_graph.size();i++){
        std::cout<<i+1<<":";
        for(int j_ = 0; j_<=i;j_++){
            std::cout<<"  ";
        }
        for(int j = i+1; j<m_graph.size();j++) {
            std::cout<<m_graph[i][j]<<",";
        }
        std::cout<<std::endl;
    }
}
/*----------------------------------------------------------------------------*/
void TopologyGraph::writeGML(){

    std::ofstream output ("TopologyGraph.gml", std::ios::out);

    output << "graph [\n";
    output << "comment \"This is a simple Topology graph with node=surface and edges between node represents connectivity of two surfaces with edge\"\n";
    output << "directed 0\n";
    output << "id 1\n";
    output << "label \"Topology graph\"\n";

    //For now the size of m_graph is the number of surfaces in the geometry
    //So we write each node of m_graph in GML format with id "i" as they are ordered from 1 to m_graph size
    for(int i = 0; i<m_graph.size();i++){
        output << "  node [\n";
        output << "    id "<<i<<"\n";
        output << "    value 1\n";
        //output << "    label \"surface "<<i<<"\"\n";
        output << "  ]\n";
    }
    //For each line of our matrix, i.e. each node i of the graph, we create an edge from node i to node j if
    //if the value of the matrix at (i,j) is not 0
    for(int i =0; i<m_graph.size();i++){
        for(int j =i; j<m_graph.size();j++){
            if(m_graph[i][j] != 0) {
                output << "  edge [\n";
                //We use add +1 to i and j because the surfaces ids start at 1 and not 0
                output << "    source " << i << "\n";
                output << "    target " << j << "\n";
                //output << "    label \"edge from " << i+1 << " to " << j+1 << "\"\n";
                //We give the value 2 if the edge is closed
                if(m_graph[i][j] == 2){
                    output << "    value 2\n";
                } else{
                    output << "    value 1\n";
                }
                output << "  ]\n";
            }
        }
    }
    output << "]";
    output<<"\n\n";

    output.close();

}
/*----------------------------------------------------------------------------*/
void TopologyGraph::colorSurfaces(){
    for(int i = 1; i <=m_graph.size(); i++){
        cad::FACSurface* surf = m_manager->getFACSurface(i);
        std::vector<Face> faces;
        surf->getMeshFaces(faces);
        //std::cout<<"FACSurface "<<i<<std::endl;
        for (auto f: faces) {
           // std::cout<<"Face id = "<<f.id()<<std::endl;
            (*m_surface_color)[f.id()] = i;
        }
    }
}
/*----------------------------------------------------------------------------*/
bool TopologyGraph::findIso(){
    //Petit hack pour tester si on peut récuperer une face selon certains critères
    //On veut une face qui n'est pas un quadrangle donc qui posède plus que 4 voisins

    for (int i = 0; i < m_graph.size(); i++) {
        //On va simplement compter le nombre d'arêtes qui bordent la face
        int cpt = 0;
        for (int j = 0; j < m_graph.size(); j++) {
            if(m_graph[i][j] == 2){
                cpt++;
            }
        }
        if(cpt == 2){
            std::cout<<"Subgraph found"<<std::endl;
            gmds::Variable<int>* BND_COLOR = m_mesh->getVariable<int,GMDS_FACE>("BND_SURFACE_COLOR");
            for(auto f: m_mesh->faces()){
                if((*BND_COLOR)[f] == i+1){
                    (*m_isomorph)[f] = 1;
                }
            }
        }
    }

}


















